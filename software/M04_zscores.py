#! /usr/bin/env python
__author__ = 'heroico'
import metax
__version__ = metax.__version__

import logging
import os
import scipy.stats as stats
import metax.KeyedDataSet as KeyedDataSet
import metax.WeightDBUtilities as WeightDBUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.Person as Person
import metax.MatrixUtilities as MatrixUtilities
import metax.ZScoreCalculation as ZScoreCalculation
import metax.Normalization as Normalization
import metax.MethodGuessing as MethodGuessing
import metax.Exceptions as Exceptions

class MResult(object):
    def __init__(self):
        self.gene = None,
        self.gene_name = None
        self.zscore = None
        self.effect_size = None
        self.p = None
        self.VAR_g = None
        self.n = None
        self.n_cov = None
        self.n_model = None
        self.gene_R2 = None
        self.gene_p = None

    HEADER="gene,gene_name,zscore,effect_size,pvalue,VAR_g,pred_perf_R2,pred_perf_p,n_snps_used,n_snps_in_cov,n_snps_in_model\n"

    def toCSVLine(self):
        line = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
               (self.gene, self.gene_name, self.zscore, self.effect_size, self.p, self.VAR_g, self.gene_R2, self.gene_p, self.n, self.n_cov, self.n_model)
        return line

class CalculateZScores(object):
    def __init__(self, args):
        self.weight_db_path = args.weight_db_path
        self.covariance = args.covariance
        self.folder_beta = args.beta_folder
        self.output_file = args.output_file
        self.selected_dosage_folder = args.selected_dosage_folder
        self.keep_ens_version = args.keep_ens_version

        self.zscore_scheme = args.zscore_scheme
        self.normalization_scheme = args.normalization_scheme

    def run(self):
        folder = os.path.split(self.output_file)[0]
        if len(folder) and not os.path.exists(folder):
            os.makedirs(folder)

        if os.path.exists(self.output_file):
            logging.info("Results path %s already exists, delete it if you want it to be calculated again", self.output_file)
            return

        people_by_id = None
        if os.path.exists(self.selected_dosage_folder):
            logging.info("Loading people")
            samples_path = Utilities.samplesInputPath(self.selected_dosage_folder)
            if samples_path is not None:
                people = Person.Person.loadPeople(samples_path)
                people_by_id = {p.id:p for p in people}

        logging.info("Loading weights from database: %s" % (self.weight_db_path))
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.weight_db_path)

        #Normalization is ignored at the moment. Not sure if it will return.
        results = None
        normalization = None
        results, normalization = self.resultsFromCovarianceFile(weight_db_logic)

        self.saveEntries(self.output_file, results)

    def resultsFromCovarianceFile(self, weight_db_logic):
        results = {}

        logging.info("Loading covariance file")
        covariance_contents = MatrixUtilities.loadMatrixFromFile(self.covariance)

        beta_contents = Utilities.contentsWithPatternsFromFolder(self.folder_beta, [])
        zscore_calculation, normalization = self.selectMethod(self.folder_beta, beta_contents, covariance_contents, weight_db_logic)

        total_entries = len(covariance_contents)
        reporter = Utilities.PercentReporter(logging.INFO, total_entries)
        i=0
        for beta_name in beta_contents:
            logging.info("Processing %s", beta_name)

            beta_path = os.path.join(self.folder_beta, beta_name)

            beta_sets = KeyedDataSet.KeyedDataSetFileUtilities.loadDataSetsFromCompressedFile(beta_path, header="")
            beta_sets = {set.name:set for set in beta_sets }
            key, check = beta_sets.iteritems().next()
            normalization.update(beta_sets)

            for gene, entry in covariance_contents.iteritems():
                weights = weight_db_logic.weights_by_gene[gene]
                process = False
                for rsid, weight in weights.iteritems():
                    if rsid in check.values_by_key:
                        process = True
                        break

                if not process:
                    logging.log(5, "No rsid in beta file for %s", gene)
                    continue

                if gene in results:
                    logging.info("Gene %s already processed", gene)
                    continue

                reporter.update(i, "%d %% of model's snp information found so far in the gwas study") #proxied by percenteage of genes

                covariance_matrix = entry[0]
                valid_rsids = entry[1]

                logging.log(7, "Calculating z score for %s", gene)

                pre_zscore, n, VAR_g, effect_size = zscore_calculation(gene, weights, beta_sets, covariance_matrix, valid_rsids)
                results[gene] = self.buildEntry(gene, weight_db_logic, weights, pre_zscore, n, VAR_g, effect_size)
                i+=1

        #second pass, for genes not in any beta file
        self.fillBlanks(results, covariance_contents, weight_db_logic, zscore_calculation)
        normalization_constant = normalization.calculateNormalization()
        return results, normalization_constant

    def buildEntry(self, gene, weight_db_logic, weights, zscore, n, VAR_g, effect_size):
        gene_entry = weight_db_logic.gene_data_for_gene[gene]
        p = "NA"
        if zscore != "NA":
            z = float(zscore)
            p = stats.norm.sf(abs(z))*2
            p = str(p)
        e = MResult()

        if not self.keep_ens_version:
            if "ENSG00" in gene and "." in gene:
                gene = gene.split(".")[0]
        e.gene = gene
        e.gene_name = gene_entry.gene_name
        e.zscore = zscore
        e.effect_size = effect_size
        e.p = p
        e.VAR_g = VAR_g
        e.n = n
        e.n_cov = str(len(weights.keys()))
        e.n_model = gene_entry.n_snp
        e.gene_R2 = gene_entry.R2
        e.gene_p = gene_entry.pval
        return e

    def fillBlanks(self, results, entries, weight_db_logic, zscore_calculation):
        dummy = { "beta_z": KeyedDataSet.KeyedDataSet(name="beta_z"),
                  "beta":KeyedDataSet.KeyedDataSet("beta"),
                  "se":KeyedDataSet.KeyedDataSet("se")  }
        for gene, entry in entries.iteritems():
            if gene in results:
                continue
            weights = weight_db_logic.weights_by_gene[gene]
            covariance_matrix = entry[0]
            valid_rsids = entry[1]

            z_score, n, VAR_g, effect_size = zscore_calculation(gene, weights, dummy, covariance_matrix, valid_rsids)
            results[gene] = self.buildEntry(gene, weight_db_logic, weights, z_score, n, VAR_g, effect_size)

    def saveEntries(self, result_path, results):
        keys = sorted(results, key=lambda key: -abs(float(results[key].zscore)) if results[key].zscore != "NA" else 0)
        with open(result_path, 'w') as file:
            file.write(MResult.HEADER)
            for key in keys:
                e = results[key]
                line = e.toCSVLine()
                file.write(line)

    def selectMethod(self, folder, beta_contents, covariance_entries, weight_db_logic):
        normalization = None
        zscore_calculation = None
        if self.zscore_scheme:
            zscore_calculation = ZScoreCalculation.ZScoreScheme(self.zscore_scheme)
            if not self.normalization_scheme:
                raise Exception("Normalization scheme is required")
        else:
            zscore_calculation, normalization = MethodGuessing.chooseZscoreSchemeFromFiles(folder, beta_contents, covariance_entries, weight_db_logic)

        if self.normalization_scheme:
            normalization = Normalization.normalizationScheme(self.normalization_scheme, covariance_entries, weight_db_logic)

        return zscore_calculation, normalization

def run(args):
    "Wrapper for common behavior for execution. "
    work = CalculateZScores(args)
    work.run()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M04_zscores.py %s: Build ZScores from GWAS data.' % (__version__,))

    parser.add_argument("--weight_db_path",
                        help="name of weight db in data folder",
                        default="data/DGN-WB_0.5.db")

    parser.add_argument("--selected_dosage_folder",
                        help="name of folder containing the selected samples, optional",
                        default="intermediate/filtered_1000GP_Phase3")

    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default="intermediate/cov/covariance.txt.gz")

    parser.add_argument("--beta_folder",
                        help="name of folder containing beta data",
                        default="intermediate/beta")

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")

    parser.add_argument("--zscore_scheme",
                        help="Scheme for zscore calculation. Options are:"
                             "'beta_z' (uses zscore of beta and sigma_l from input file), default;"
                            " 'beta_z_and_ref' (uses zscore of beta and sigma_l from reference population);"
                            " 'metaxcan' ('bare metal' MetaXcan, normalization recommended); "
                            " 'metaxcan_from_reference' ('bare metal' MetaXcan, using reference variance);",
                        default=None)

    parser.add_argument("--normalization_scheme",
                        help="Scheme for zscore normalization, relevant for 'global normalization' scheme. Options are:"
                            "'none';"
                            " 'from_pheno', estimate normalization constant from phenotype file, needs 'sigma_l' and 'standard error' in phenotype;"
                            " 'from_reference', estimate normalization constant from reference, needs 'standard error' on phenotype",
                        default=None)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--keep_ens_version",
                        help="If set, will keep the -version- postfix in gene id.",
                    action="store_true",
                    default=False)

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)


    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))
    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except NameError as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)
        except Exceptions.ReportableException, e:
            logging.error(e.msg)

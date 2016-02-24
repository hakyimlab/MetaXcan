#! /usr/bin/env python
__author__ = 'heroico'

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


class CalculateZScores(object):
    def __init__(self, args):
        self.weight_db_path = args.weight_db_path
        self.covariance = args.covariance
        self.folder_beta = args.beta_folder
        self.output_file = args.output_file
        self.selected_dosage_folder = args.selected_dosage_folder

        self.zscore_scheme = args.zscore_scheme
        self.normalization_scheme = args.normalization_scheme

    def run(self):
        folder = os.path.split(self.output_file)[0]
        if len(folder) and not os.path.exists(folder):
            os.makedirs(folder)

        if os.path.exists(self.output_file):
            logging.info("Results path %s already exists, delete it if you want it to be calculated again", self.output_file)
            return

        logging.info("Loading people")
        people_by_id = None
        if os.path.exists(self.selected_dosage_folder):
            samples_path = Utilities.samplesInputPath(self.selected_dosage_folder)
            if samples_path is not None:
                people_by_id = Person.Person.peopleByIdFromFile(samples_path)

        logging.info("Loading weight db")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.weight_db_path)

        results = None
        normalization = None
        results, normalization = self.resultsFromCovarianceFile(weight_db_logic)

        self.saveEntries(self.output_file, results, normalization)

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
                weights = weight_db_logic.weights_by_gene_name[gene]
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

                reporter.update(i, "%d done")

                covariance_matrix = entry[0]
                valid_rsids = entry[1]

                logging.log(7, "Calculating z score for %s", gene)
                pre_zscore, n, VAR_g = zscore_calculation(gene, weights, beta_sets, covariance_matrix, valid_rsids)
                results[gene] = self.buildEntry(gene, weight_db_logic, weights, pre_zscore, n, VAR_g)
                i+=1

        #second pass, for genes not in any beta file
        self.fillBlanks(results, covariance_contents, weight_db_logic, zscore_calculation)
        normalization_constant = normalization.calculateNormalization()
        return results, normalization_constant

    def buildEntry(self, gene, weight_db_logic, weights, zscore, n, VAR_g):
        gene_entry = weight_db_logic.gene_data_for_gene[gene]
        p = "NA"
        if zscore != "NA":
            z = float(zscore)
            p = stats.norm.sf(abs(z))*2
            p = str(p)
        return (gene, gene_entry.gene_name, zscore, VAR_g, n, str(len(weights.keys())), gene_entry.n_snp, gene_entry.R2, p, )

    def fillBlanks(self, results, entries, weight_db_logic, zscore_calculation):
        dummy = { "beta_z": KeyedDataSet.KeyedDataSet(name="beta_z"),
                  "beta":KeyedDataSet.KeyedDataSet("beta"),
                  "se":KeyedDataSet.KeyedDataSet("se")  }
        for gene, entry in entries.iteritems():
            if gene in results:
                continue
            weights = weight_db_logic.weights_by_gene_name[gene]
            covariance_matrix = entry[0]
            valid_rsids = entry[1]

            z_score, n, VAR_g = zscore_calculation(gene, weights, dummy, covariance_matrix, valid_rsids)
            results[gene] = self.buildEntry(gene, weight_db_logic, weights, z_score, n, VAR_g)

    def saveEntries(self, result_path, results, normalization):
        keys = sorted(results, key=lambda key: -abs(float(results[key][2])) if results[key][2] != "NA" else 0)
        with open(result_path, 'w') as file:
            file.write("gene,gene_name,zscore,pvalue,pred_perf_R2,VAR_g,n,covariance_n,model_n\n")
            for key in keys:
                entry = results[key]
                gene = entry[0]
                gene_name = entry[1]
                z_score = entry[2]
                if normalization != 1:
                    if z_score != "NA":
                        z = float(z_score)
                        z_score = str(z*normalization)
                VAR_g = entry[3]
                n = entry[4]
                covariance_n = entry[5]
                model_n = entry[6]
                R2 = entry[7]
                p = entry[8]
                line = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (gene, gene_name, z_score, p, R2, VAR_g, n, covariance_n, model_n)
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
    if args.throw:
        work.run()
    else:
        try:
            work.run()
        except Exception as e:
            logging.info("Unexpected error: %s", str(e))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build betas from GWAS data.')

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

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))
    run(args)



#! /usr/bin/env python
__author__ = 'heroico'

import logging
import os
import re
import scipy.stats as stats
import metax.KeyedDataSet as KeyedDataSet
import metax.WeightDBUtilities as WeightDBUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.DBLoaders as DBLoaders
import metax.Person as Person
import metax.Formats as Formats
import metax.MatrixUtilities as MatrixUtilities
import metax.ZScoreCalculation as ZScoreCalculation
import metax.Normalization as Normalization
import metax.MethodGuessing as MethodGuessing


class CalculateZScores(object):
    def __init__(self, args):
        self.weight_db_path = args.weight_db_path
        self.folder_covariance = args.covariance_folder
        self.folder_beta = args.beta_folder
        self.output_file = args.output_file
        self.selected_dosage_folder = args.selected_dosage_folder

        self.zscore_scheme = args.zscore_scheme
        self.normalization_scheme = args.normalization_scheme

        self.input_format = args.input_format

        self.gene_in_covariance_regexp = re.compile(args.covariance_file_pattern)

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
        if self.input_format == Formats.DatabaseFile:
            results, normalization = self.resultsFromDatabases(weight_db_logic)
        elif self.input_format == Formats.FlatFile:
            results, normalization = self.resultsFromCovarianceFile(weight_db_logic)
        else:
            logging.info("Wrong input format: %s", self.input_format)
            assert False

        self.saveEntries(self.output_file, results, normalization)

    def resultsFromCovarianceFile(self, weight_db_logic):
        results = {}

        logging.info("Loading covariance file")
        file = Utilities.contentsWithPatternsFromFolder(self.folder_covariance, [".gz"])[0]
        path = os.path.join(self.folder_covariance, file)
        covariance_contents = MatrixUtilities.loadMatrixFromFile(path)

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

    def resultsFromDatabases(self, weight_db_logic):
        results = {}

        covariance_contents = Utilities.contentsWithPatternsFromFolder(self.folder_covariance, ["cov"])
        beta_contents = Utilities.contentsWithPatternsFromFolder(self.folder_beta, ["beta"])
        zscore_calculation, normalization = self.selectMethod(self.folder_beta, beta_contents, covariance_contents, weight_db_logic)

        CHROMOSOMES = ["chr"+str(x) for x in xrange(1, 23)]
        for chromosome in CHROMOSOMES:
            logging.info("Processing %s", chromosome)
            chromosome_covariances = Utilities.removeNamesWithPatterns(covariance_contents, [chromosome+"-"])

            beta_name = Utilities.removeNameWithPatterns(beta_contents, [chromosome+"."])
            if beta_name:
                logging.info("beta_name %s", beta_name)
            beta_path = os.path.join(self.folder_beta, beta_name)
            beta_sets = KeyedDataSet.KeyedDataSetFileUtilities.loadDataSetsFromCompressedFile(beta_path, header="")
            normalization.update(beta_sets)

            total_covariances = len(chromosome_covariances)
            reporter = Utilities.PercentReporter(logging.INFO, total_covariances)
            for i,chromosome_covariance in enumerate(chromosome_covariances):
                reporter.update(i, "Chromosome %s, %s" %(chromosome_covariance, "%d done"))

                regexp_search = self.gene_in_covariance_regexp.search(chromosome_covariance)
                gene = regexp_search.group(1)
                weights = weight_db_logic.weights_by_gene_name[gene]
                if gene in results:
                    logging.info("Gene %s already done, skipping")
                    continue

                logging.log(7, "Loading covariance for %s", gene)
                covariance_matrix, valid_rsids = self.loadCovarianceFromDB(weights, chromosome_covariance)

                logging.log(7, "Calculating z score for %s", gene)
                pre_zscore, n, VAR_g = zscore_calculation(gene, weights, beta_sets, covariance_matrix, valid_rsids)
                results[gene] = self.buildEntry(gene, weight_db_logic, weights, pre_zscore, n, VAR_g)

        normalization_constant = normalization.calculateNormalization()
        return results, normalization_constant

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

    def loadCovarianceFromDB(self, weights, chromosome_covariance_name):
        rsids = weights.keys()
        db_path = os.path.join(self.folder_covariance, chromosome_covariance_name)
        covariance_matrix, valid_rsids = DBLoaders.DBLoaders.loadCovarianceMatrix(db_path, rsids)
        return covariance_matrix, valid_rsids

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

    parser.add_argument("--covariance_folder",
                        help="name of folder containing covariance data",
                        default="intermediate/cov")

    parser.add_argument("--covariance_file_pattern",
                        help="pattern for covariance database file names, regexp to select name of gene. Optional, not used if flat file format.",
                        default='cov-1000GP_Phase3_chr(?<!\d)\d{1,2}-(.*)-DGN-WB_0.5.db')

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

    parser.add_argument('--input_format',
                   help='Input covariance files format. Valid options are: CovarianceDatabase, CovarianceFile',
                   default=Formats.FlatFile)

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



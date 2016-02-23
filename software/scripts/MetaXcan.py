#! /usr/bin/env python

__author__ = 'heroico'

__version__ = 0.1
from subprocess import call
import logging
import metax.Logging as Logging
import metax.ZScoreCalculation as ZScoreCalculation
import metax.Formats as Formats

class MetaXcanProcess(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        rcode = self.buildBetas()
        if rcode != 0:
            logging.info("Beta step failed!")
            return

        self.buildZScores()

    def buildBetas(self):
        logging.info("Processing betas!")
        command = "./M03_betas.py"
        command += " --weight_db_path " + self.args.weight_db_path
        command += " --gwas_folder " + self.args.gwas_folder
        command += " --output_folder " + self.args.beta_folder
        command += " --verbosity " + self.args.verbosity
        if self.args.gwas_file_pattern:
            command += " --gwas_file_pattern " + self.args.gwas_file_pattern
        command += " --a1_column " + self.args.a1_column
        command += " --a2_column " + self.args.a2_column
        command += " --snp_column " + self.args.snp_column
        if self.args.scheme:
            command += " --scheme " + self.args.scheme
        if self.args.or_column:
            command += " --or_column " + self.args.or_column
        if self.args.beta_column:
            command += " --beta_column " + self.args.beta_column
        if self.args.beta_sign_column:
            command += " --beta_sign_column " + self.args.beta_sign_column
        if self.args.beta_zscore_column:
            command += " --beta_zscore_column " + self.args.beta_zscore_column
        if self.args.se_column:
            command += " --se_column " + self.args.se_column
        if self.args.frequency_column:
            command += " --frequency_column " + self.args.frequency_column
        if self.args.pvalue_column:
            command += " --pvalue_column " + self.args.pvalue_column
        if self.args.compressed:
            command += " --compressed"
        if self.args.throw:
            command += " --throw"
        if self.args.separator:
            command = " --separator "+self.args.separator
        rcode = call(command.split(" "))
        return rcode

    def buildZScores(self):
        logging.info("Calculating ZScores!")
        command = "./M04_zscores.py"
        command += " --weight_db_path " + self.args.weight_db_path
        command += " --selected_dosage_folder " + self.args.selected_dosage_folder
        command += " --covariance_folder " + self.args.covariance_folder
        command += " --covariance_file_pattern " + self.args.covariance_file_pattern
        command += " --beta_folder " + self.args.beta_folder
        command += " --output_file " + self.args.output_file
        command += " --input_format " + self.args.covariance_input_format
        command += " --verbosity " + self.args.verbosity
        if self.args.zscore_scheme:
            command += " --zscore_scheme " + self.args.zscore_scheme
        if self.args.normalization_scheme:
            command += " --normalization_scheme " + self.args.normalization_scheme
        if self.args.throw:
            command += " --throw"
        rcode = call(command.split(" "))
        return rcode

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Will estimate MetaXcan results from a set of snp covariance matrices, a model database, and GWAS beta files.')

#weight db model
    parser.add_argument("--weight_db_path",
                        help="name of weight db in data folder",
                        default="data/DGN-WB_0.5.db")

#GWAS betas
    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    parser.add_argument("--scheme",
                        help="Type of beta data preprocessing, optional. Options are: "
                        "'beta' (provide (beta or OR));"
                        "'beta_se' (provide (beta or OR) and standard error); "
                        "'beta_se_to_z' (provide (beta or OR) and standard error), and Zscore of beta will be output;"
                        "'z' (provide zscore of beta),"
                        " 'beta_sign_p' (sign of beta, and pvalue); beta_p (beta and pvalue)",
                        default=None)

    parser.add_argument("--or_column",
                    help="Name of column containing Odd Ratios in input files. Either 'OR_column' or 'beta_column' must be provided",
                    default=None)

    parser.add_argument("--pvalue_column",
                    help="Name of column containing p-value in input files.",
                    default=None)

    parser.add_argument("--beta_sign_column",
                    help="Name of column containing sign of beta in input files.",
                    default=None)

    parser.add_argument("--beta_column",
                    help="Name of column containing betas in input files. Either 'OR_column' or 'beta_column' must be provided",
                    default=None)

    parser.add_argument("--se_column",
                    help="Name of column containing standard error in input file",
                    default=None)

    parser.add_argument("--beta_zscore_column",
                    help="Name of column containing beta's zscore in input file.",
                    default=None)

    parser.add_argument("--frequency_column",
                    help="Name of column containing frequency in input file",
                    default="FRQ")

    parser.add_argument("--a1_column",
                    help="Name of column containing allele 1 in input file",
                    default="A1")

    parser.add_argument("--a2_column",
                    help="Name of column containing allele 2 in input file",
                    default="A2")

    parser.add_argument("--snp_column",
                    help="Name of column containing snp in input file",
                    default="SNP")

    parser.add_argument("--compressed",
                    help="Wether input files are gzip compressed file",
                    action="store_true",
                    default=False)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)

#both

    parser.add_argument("--beta_folder",
                        help="name of folder to put beta parsing results in",
                        #default="intermediate/beta")
                        default="intermediate/beta")

# ZScore calculation
    parser.add_argument("--selected_dosage_folder",
                        help="name of folder containing the selected samples, optional",
                        #default="intermediate/filtered_1000GP_Phase3")
                        default="intermediate/filtered_1000GP_Phase3")

    parser.add_argument("--covariance_folder",
                        help="name of folder containing covariance data",
                        #default="intermediate/1000GP_Phase3_chr_cov")
                        default="intermediate/cov")

    parser.add_argument("--covariance_file_pattern",
                        help="pattern for covariance file names, must contain 'EXTRACT_GENE'",
                        default='cov-1000GP_Phase3_chr(?<!\d)\d{1,2}-(.*)-DGN-WB_0.5.db')

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

    parser.add_argument('--covariance_input_format',
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

    work = MetaXcanProcess(args)
    work.run()
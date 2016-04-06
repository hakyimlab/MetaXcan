#! /usr/bin/env python

__author__ = 'heroico'

import metax
__version__ = metax.__version__
import logging
import metax.Logging as Logging
import M03_betas
import M04_zscores
import metax.Exceptions as Exceptions

class MetaXcanProcess(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        self.buildBetas()
        self.buildZScores()

    def buildBetas(self):
        logging.info("Processing betas!")
        self.args.output_folder = args.beta_folder
        M03_betas.run(self.args)

    def buildZScores(self):
        logging.info("Calculating ZScores!")
        M04_zscores.run(self.args)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MetaXcan.py %s:  Will estimate MetaXcan results from a set of snp covariance matrices, a model database, and GWAS beta files.' % (__version__))

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
                    default=None)

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

    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.",
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

    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        #default="intermediate/1000GP_Phase3_chr_cov")
                        default="intermediate/cov/covariance.txt.gz")

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

    parser.add_argument("--keep_ens_version",
                    help="If set, will keep the -version- postfix in gene id.",
                    action="store_true",
                    default=False)


    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    try:
        work = MetaXcanProcess(args)
        work.run()
    except Exceptions.ReportableException, e:
        logging.error(e.msg)

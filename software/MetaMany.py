#! /usr/bin/env python

__author__ = 'heroico, Eric Torstenson'

"""This modified version of MetaXcan provides support for analyzing multiple
   tissues at with a single run, which is purely for convenience.

   Big differences include users specifying one or more tissue databases,
   and directing the script to the covariance database directory, working
   under the assumption that the covariances will be named nearly identically
   to the tissue database, with the exception of extension (and path). This
   is true with the current build of GTEx databases found on the box
   directory.

   Additional parameter changes include:
    * Specifying multiple weight databases
    * Specifying covariance directory, rather than individual files
    * Specifying output directory, rather than named file (results will be
      named identically to the database files, except for extension and
      path.
    * ? (probably others)
   """
import metax
__version__ = metax.__version__
import logging
import metax.Logging as Logging
import M03_betas
import M04_zscores
import metax.Exceptions as Exceptions
import os
import metax.WeightDBUtilities as WeightDBUtilities
import metax.Utilities as Utilities

class MetaXcanProcess(object):
    def __init__(self, args):
        self.args = args


    def buildBetas(self, db_filename):
        filebase = os.path.basename(filename).replace(".db", "")
        output_folder = self.args.output_folder
        logging.info("Processing betas for %s" % (db_filename))
        self.args.weight_db_path = os.path.abspath(db_filename)
        self.args.covariance = os.path.join(self.args.covariance_directory, filebase) + ".txt.gz"
        self.args.output_file = os.path.join(self.args.output_directory, filename) + ".csv"

        logging.info("Loading weight model")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.args.weight_db_path)

        names = Utilities.contentsWithRegexpFromFolder(self.args.gwas_folder, self.args.gwas_regexp)

        if not os.path.exists(self.args.output_folder):
            os.makedirs(self.args.output_folder)

        for name in names:
            try:
                resultName = "%s-%s" % (name.split("/")[-1].split(".")[0], filebase)
                self.buildBetas(weight_db_logic,name,resultName)
            # This just means that there is some extra stuff inside that directory,
            # so I'm thinking we want to ignore it.
            except Exceptions.BadFilename as e:
                logging.info("Wrong file name: %s, skipping", e.msg)
                pass
    def run(self):
        self.args.output_folder = args.beta_folder
        for weight_db in self.args.weight_dbs:
            filename = weight_db.name
            filebase = os.path.basename(filename).replace(".db", "")
            weight_db.close()
            logging.info("Processing betas for %s" % (filebase))
            self.args.weight_db_path = os.path.abspath(filename)
            self.args.covariance = os.path.join(self.args.covariance_directory, filebase) + ".txt.gz"

            self.args.output_file = os.path.join(self.args.output_directory, filename) + ".csv"
            M03_betas.run(self.args)

            logging.info("Calculating ZScores for %s" % (filebase))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MetaXcan.py %s:  Will estimate MetaXcan results from a set of snp covariance matrices, a model database, and GWAS beta files.' % (__version__))

#weight db model
    parser.add_argument('weight_dbs', metavar='DB', type=argparse.FileType('r'),
                        nargs='+', help="weight dbs to be used")

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

    parser.add_argument("--covariance_directory",
                        help="directory where covariance files can be found",
                        default="intermediate/cov")

    parser.add_argument("--output_directory",
                        help="name of output file",
                        default="results")

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

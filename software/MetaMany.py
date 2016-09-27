#! /usr/bin/env python

import metax
import sys
__version__ = "1.0x" + metax.__version__
import logging
import metax.Logging as Logging
import M03_betas
import M04_zscores
import metax.Exceptions as Exceptions
import os
import metax.WeightDBUtilities as WeightDBUtilities
import metax.Utilities as Utilities
import re
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

   Changes in 1.0.0.3
    * Added a version prefix to disambiguate changes in the script to the
      underlying library.
    * To provide more flexibility with tissue database structure, the command
      line interface has changed slightly. For the August 2016 release,
      the default values should work fine. So, you won't need to provide
      --covariance_directory (or if you want to explicitly provide defaults
      for documentation in your scripts, set it to be SAME instead of a
      PATH name)
    * Betas are now written into a subdirectory inside the output directory
      (these are written inside directories named after the tissues) Results
      are written directly inside the output directory and are prefixed as
      follows:
        GWAS-Tissue.csv[.gz]
      Where GWAS is the first gwas filename found inside the gwas directory
      Tissue is the tissue for which the zscores are calculated
      and the optional .gz is there if the --compressed option is set
   """


class MetaXcanProcess(object):
    def __init__(self, args):
        self.args = args
        self.gwas_regexp = None
        if args.gwas_file_pattern:
            self.gwas_regexp = re.compile(args.gwas_file_pattern)

    def buildBetas(self, db_filename):
        filebase = os.path.basename(db_filename).replace(".db", "")
        output_folder = os.path.abspath(self.args.output_directory)

        logging.info("Processing betas for %s" % (db_filename))
        self.args.weight_db_path = os.path.abspath(db_filename)
        cov_directory = self.args.covariance_directory
        if cov_directory.upper() == "SAME":
            cov_directory = "/".join(self.args.weight_db_path.split("/")[0:-1])

        extComponents = self.args.covariance_suffix.split("..")

        if len(extComponents) > 1:
            covext = "..".join(extComponents[0:-1])
            dbext = extComponents[-1]
            filebase = db_filename.replace(dbext, "")
            self.args.covariance = "%s/%s%s" % (cov_directory, filebase.split("/")[-1], covext)
        else:
            self.args.covariance = "%s/%s%s" % (
            cov_directory, filebase.strip("/")[-1], self.args.covariance_suffix)
        file_prefix = filebase.split("/")[-1].split(".")[0]
        beta_output = os.path.join(output_folder, file_prefix)
        logging.info("Writing betas to %s" % (beta_output))

        self.args.output_folder = beta_output

        logging.info("Loading weight model")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.args.weight_db_path)

        betaScript = M03_betas.GetBetas(self.args)
        names = Utilities.contentsWithRegexpFromFolder(self.args.gwas_folder, betaScript.gwas_regexp)

        if not os.path.exists(beta_output):
            os.makedirs(beta_output)
        betaScript.output_folder = beta_output              #os.path.join(output_folder, filebase)
        if not os.path.exists(betaScript.output_folder):
            os.makedirs(betaScript.output_folder)

        report_prefix = None
        for name in names:
            if report_prefix is None:
                report_prefix = name.split("/")[-1].split(".")[0]
            try:
                betaScript.buildBetas(weight_db_logic,name)

            # This just means that there is some extra stuff inside that directory,
            # so I'm thinking we want to ignore it.
            except Exceptions.BadFilename as e:
                logging.info("Wrong file name: %s, skipping", e.msg)
                pass

        suffix = ".csv"
        if args.compressed:
            suffix += ".gz"
        self.args.output_file = os.path.join(output_folder,
                                             report_prefix + "-" + file_prefix + suffix)  # output_folder       #os.path.join(output_folder, file_prefix) + ".csv"

        # ZScores
        logging.info("Calculating ZScores for %s" % (filebase))
        zscoreScript = M04_zscores.CalculateZScores(self.args)
        zscoreScript.folder_beta = betaScript.output_folder
        zscoreScript.run()

    def run(self):
        self.args.output_folder = args.beta_folder
        for weight_db in self.args.weight_dbs:
            weight_db.close()
            self.buildBetas(weight_db.name)




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MetaMany.py %s:  Automates the execution of MetaXcan over multiple tissues. ' % (__version__),
                                     epilog="""
MetaMany makes certain assumptions about the path relationships between the
tissue databases and their corresponding covariance files. While the defaults
are intended to relate to the current databases' organization and naming, it
is possible that users will have to make changes if either the database file
hierarchy is different (and newer) than this script, or when trying to use
older versions (or perhaps those created by the end user which don't follow the
same organization. In these cases, users should pay careful attention to the
arguments --covariance_directory and --covariance_suffix. """)

    parser.add_argument("-v", "--version", help="Report the version",
                        action="store_true",
                        default=False)
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

    parser.add_argument("--non_effect_allele_column",
                    help="Name of column containing non-effect allele in input file ('reference allele', if following PrediXcan format, and plink --dosage format philosophy)",
                    default="A2")

    parser.add_argument("--effect_allele_column",
                    help="Name of column containing effect (or dosage) allele in input file (dosage/effect allele)",
                    default="A1")

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
    # Added to support GWAS utilities
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

    parser.add_argument("--covariance_directory",
                        help="directory where covariance files can be found (or SAME if covariance sits beside the .db file",
                        default="SAME")

    parser.add_argument("--covariance_suffix",
                        help="Suffix associated with the covariate files. covext-dbext (where ..dbext is the portion of the db file to be replaced by the coviarance extention. )",
                        default=".txt.gz.._0.5.db")

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

    if "-v" in sys.argv or "--verbose" in sys.argv:
        print metax.__version__
        sys.exit(0)
    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    try:
        work = MetaXcanProcess(args)
        work.run()
    except Exceptions.ReportableException, e:
        logging.error(e.msg)

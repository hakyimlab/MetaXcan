#! /usr/bin/env python

import metax
import sys
__version__ = "1.0x" + metax.__version__
import logging
import re
import os
from metax import Logging
from metax import Exceptions
from metax import Utilities
from metax.gwas import Utilities as GWASUtilities

import SPrediXcan

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

   Changes in 1.0.0.5:
    * Refactor, cleaned up stale code, and dropped deprecated parameters

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
      and the optional .gz is there if the --compressed_gwas option is set
   """

def get_name_prefix(args):
    if args.gwas_folder:
        regexp = re.compile(args.gwas_file_pattern) if args.gwas_file_pattern else None
        names = Utilities.contentsWithRegexpFromFolder(args.gwas_folder, regexp)
        name = names[0]
        report_prefix = get_result_prefix(args, name)
    else:
        report_prefix = get_result_prefix(args, args.gwas_file)
    return report_prefix

def get_result_prefix(args, name):
    if args.output_file_prefix:
        return args.output_file_prefix
    n = os.path.splitext(args.gwas_file)[0] if "." in args.gwas_file else args.gwas_file
    return n

def process(args, db_filename):
    filebase = os.path.basename(db_filename).replace(".db", "")
    output_folder = os.path.abspath(args.output_directory)

    logging.info("Processing %s" % (db_filename))
    args.model_db_path = os.path.abspath(db_filename)
    cov_directory = args.covariance_directory
    if cov_directory.upper() == "SAME":
        cov_directory = "/".join(args.model_db_path.split("/")[0:-1])

    extComponents = args.covariance_suffix.split("..")

    if len(extComponents) > 1:
        covext = "..".join(extComponents[0:-1])
        dbext = extComponents[-1]
        filebase = db_filename.replace(dbext, "")
        args.covariance = "%s/%s%s" % (cov_directory, filebase.split("/")[-1], covext)
    else:
        args.covariance = "%s/%s%s" % (
        cov_directory, filebase.strip("/")[-1], args.covariance_suffix)
    file_prefix = filebase.split("/")[-1].split(".")[0]

    suffix = ".csv"
    report_prefix = get_name_prefix(args)
    args.output_file = os.path.join(output_folder, report_prefix + "-" + file_prefix + suffix)  # output_folder       #os.path.join(output_folder, file_prefix) + ".csv"

    # Run!
    SPrediXcan.run(args)

def run(args):
    for weight_db in args.weight_dbs:
        weight_db.close()
        process(args, weight_db.name)

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

    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")

    parser.add_argument("-v", "--version", help="Report the version", action="store_true", default=False)
#weight db model
    parser.add_argument('weight_dbs', metavar='DB', type=argparse.FileType('r'), nargs='+', help="weight dbs to be used")

#GWAS betas
    parser.add_argument("--gwas_file", help="Load a single GWAS file. (Alternative to providing a gwas_folder and gwas_file_pattern)")

    parser.add_argument("--gwas_folder", help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.")
    parser.add_argument("--gwas_file_pattern", help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).")

    GWASUtilities.add_gwas_arguments_to_parser(parser)

# ZScore calculation
    parser.add_argument("--stream_covariance", help="Option to better handle large covariances, slower but less memory consuming", action="store_true")
    parser.add_argument("--single_snp_model", action="store_true", help="Models are comprised of a single snp per gene", default=False)
    parser.add_argument("--covariance_directory", help="directory where covariance files can be found (or SAME if covariance sits beside the .db file", default="SAME")
    parser.add_argument("--covariance_suffix", help="Suffix associated with the covariance files. covext-dbext (where ..dbext is the portion of the db file to be replaced by the coviarance extension. )", default=".txt.gz.._0.5.db")
    parser.add_argument("--output_directory", help="name of output file", default="results")
    parser.add_argument("--output_file_prefix", help="name of output file", default="results")
    parser.add_argument("--additional_output", help="If set, will output additional information.", action="store_true", default=False)
    parser.add_argument("--remove_ens_version", help="If set, will keep the -version- postfix in gene id.", action="store_true", default=False)
    parser.add_argument("--overwrite", help="If set, will overwrite the results file if it exists.", action="store_true", default=False)
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--MAX_R", help="Run only for the first R genes", type=int, default=None)


    if "-v" in sys.argv or "--verbose" in sys.argv:
        print(metax.__version__)
        sys.exit(0)
    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException as e:
            logging.error(e.msg)

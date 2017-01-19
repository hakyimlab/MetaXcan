#! /usr/bin/env python

__author__ = 'heroico'

import metax
__version__ = metax.__version__
import logging
import os

from metax import Exceptions
from metax import Logging
from metax.gwas import Utilities as GWASUtilities

import M03_betas
import M04_zscores


def run(args):
    if not args.overwrite and os.path.exists(args.output_file):
        logging.info("%s already exists, move it or delete it if you want it done again", args.output_file)
        return
    if not args.model_db_path:
        logging.info("Need to provide a model database file path")
        return
    args.output_folder = None
    g = M03_betas.run(args)
    M04_zscores.run(args, g)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MetaXcan.py %s:  Will estimate MetaXcan results from a set of snp covariance matrices, a model database, and GWAS beta files.' % (__version__))

#weight db model
    parser.add_argument("--model_db_path",
                        help="name of model db in data folder",
                        default=None)

#GWAS betas
    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    GWASUtilities.add_gwas_arguments_to_parser(parser)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)

    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.",
                        default=None)

# ZScore calculation
    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        #default="intermediate/1000GP_Phase3_chr_cov")
                        default="intermediate/cov/covariance.txt.gz")

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")

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

    parser.add_argument("--overwrite",
                        help="If set, will overwrite the results file if it exists.",
                    action="store_true",
                    default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException, e:
            logging.error(e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)

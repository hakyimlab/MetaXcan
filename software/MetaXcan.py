#! /usr/bin/env python

__author__ = 'heroico'

import metax
__version__ = metax.__version__
import logging
import os
import copy

from metax import Exceptions
from metax import Logging
from metax.gwas import Utilities as GWASUtilities

import M03_betas
import M04_zscores


def run(args):
    if not args.overwrite and args.output_file and os.path.exists(args.output_file):
        logging.info("%s already exists, move it or delete it if you want it done again", args.output_file)
        return
    if not args.model_db_path:
        logging.info("Need to provide a model database file path")
        return
    M03_args = copy.copy(args)
    M03_args.output_folder = None
    M03_args.output = None
    g = M03_betas.run(M03_args)

    M04_zscores.run(args, g)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MetaXcan.py %s:  Will estimate MetaXcan results from a set of snp covariance matrices, a model database, and GWAS beta files.' % (__version__))

#weight db model
    parser.add_argument("--model_db_path", help="name of model db in data folder")
    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")
#GWAS betas
    parser.add_argument("--gwas_file", help="Load a single GWAS file. (Alternative to providing a gwas_folder and gwas_file_pattern)")

    parser.add_argument("--gwas_folder", help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.")
    parser.add_argument("--gwas_file_pattern", help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).")

    GWASUtilities.add_gwas_arguments_to_parser(parser)

# ZScore calculation
    parser.add_argument("--single_snp_model", action="store_true", help="Models are comprised of a single snp per gene", default=False)
    parser.add_argument("--covariance", help="name of file containing covariance data")
    parser.add_argument("--stream_covariance", help="Option to better handle large covariances, slower but less memory consuming", action="store_true")
    parser.add_argument("--output_file", help="name of output file")
    parser.add_argument("--remove_ens_version", help="If set, will drop the -version- postfix in gene id.", action="store_true", default=False)
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--overwrite", help="If set, will overwrite the results file if it exists.", action="store_true", default=False)
    parser.add_argument("--additional_output", help="If set, will output additional information.", action="store_true", default=False)
    parser.add_argument("--MAX_R", help="Run only for the first R genes", type=int, default=None)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException as e:
            logging.error(e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)

#! /usr/bin/env python
__author__ = 'heroico'
import metax
__version__ = metax.__version__

import logging
import os

from timeit import default_timer as timer

from metax import Logging
from metax import Utilities
from metax import Exceptions
from metax.metaxcan import AssociationCalculation
from metax.metaxcan import Utilities as MetaxcanUtilities


def run(args, _gwas=None):
    start = timer()
    if not args.overwrite and os.path.exists(args.output_file):
        logging.info("%s already exists, move it or delete it if you want it done again", args.output_file)
        return
    logging.info("Started metaxcan association")

    context = MetaxcanUtilities.build_context(args, _gwas)

    model_snps = context.get_model_snps()
    total_snps = len(model_snps)
    snps_found=set()
    reporter = Utilities.PercentReporter(logging.INFO, total_snps)

    i_genes, i_snps = context.get_data_intersection()

    results = []
    for gene in i_genes:
        r, snps = AssociationCalculation.association(gene, context, return_snps=True)
        results.append(r)
        snps_found.update(snps)
        reporter.update(len(snps_found), "%d %% of model's snps found so far in the gwas study")

    Utilities.ensure_requisite_folders(args.output_file)

    reporter.update(len(snps_found), "%d %% of model's snps used", force=True)
    results = AssociationCalculation.dataframe_from_results(zip(*results))
    results = MetaxcanUtilities.format_output(results, context, args.remove_ens_version)
    results.to_csv(args.output_file, index=False)
    end = timer()
    logging.info("Sucessfully processed metaxcan association in %s seconds"%(str(end - start)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M04_zscores.py %s: Build ZScores from GWAS data.' % (__version__,))

    parser.add_argument("--model_db_path",
                        help="name of weight db in data folder",
                        default=None)

    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default=None)

    parser.add_argument("--beta_folder",
                        help="name of folder containing beta data",
                        default=None)

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--remove_ens_version",
                        help="If set, will drop the -version- postfix in gene id.",
                    action="store_true",
                    default=False)

    parser.add_argument("--overwrite",
                        help="If set, will overwrite the results file if it exists.",
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
        except Exceptions.ReportableException as e:
            logging.error("Error:%s", e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)

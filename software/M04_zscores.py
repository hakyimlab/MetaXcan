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


def run_metaxcan(args, context):
    logging.info("Started metaxcan association")
    model_snps = context.get_model_snps()
    total_snps = len(model_snps)
    snps_found=set()
    reporter = Utilities.PercentReporter(logging.INFO, total_snps)

    i_genes, i_snps = context.get_data_intersection()

    results = []
    additional = []
    for i,gene in enumerate(i_genes):
        if args.MAX_R and i+1>args.MAX_R:
            logging.log("Early exit condition met")
            break
        logging.log(9, "Processing gene %i:%s", i, gene)
        r, snps = AssociationCalculation.association(gene, context, return_snps=True)
        results.append(r)
        snps_found.update(snps)
        reporter.update(len(snps_found), "%d %% of model's snps found so far in the gwas study")
        if args.additional_output:
            stats_ = AssociationCalculation.additional_stats(gene, context)
            additional.append(stats_)

    reporter.update(len(snps_found), "%d %% of model's snps used", force=True)

    results = AssociationCalculation.dataframe_from_results(results)
    results = MetaxcanUtilities.format_output(results, context, args.remove_ens_version)

    if args.additional_output:
        additional = AssociationCalculation.dataframe_from_aditional_stats(additional)
        results = MetaxcanUtilities.merge_additional_output(results, additional, context, args.remove_ens_version)

    if args.output_file:
        Utilities.ensure_requisite_folders(args.output_file)
        results.to_csv(args.output_file, index=False, na_rep="NA")

    return results

def run(args, _gwas=None):
    if not args.output_file and not args.additional_output:
        logging.info("Provide at least --output_file or --additional_output")

    if args.output_file and not args.overwrite and os.path.exists(args.output_file):
        logging.info("%s already exists, move it or delete it if you want it done again", args.output_file)
        return

    start = timer()
    logging.info("Started metaxcan process")

    context = MetaxcanUtilities.build_context(args, _gwas)

    results = run_metaxcan(args, context)

    end = timer()
    logging.info("Sucessfully processed metaxcan association in %s seconds"%(str(end - start)))
    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M04_zscores.py %s: Build ZScores from GWAS data.' % (__version__,))

    parser.add_argument("--model_db_path", help="name of weight db in data folder")
    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")
    parser.add_argument("--covariance", help="name of file containing covariance data")
    parser.add_argument("--stream_covariance", help="Option to better handle large covariances, slower but less memory consuming", action="store_true")
    parser.add_argument("--beta_folder", help="name of folder containing GWAS effect data")
    parser.add_argument("--output_file", help="name of output file")
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--MAX_R", help="Run only for the first R genes", type=int, default=None)
    parser.add_argument("--remove_ens_version", help="If set, will drop the -version- postfix in gene id.", action="store_true", default=False)
    parser.add_argument("--overwrite", help="If set, will overwrite the results file if it exists.", action="store_true", default=False)
    parser.add_argument("--additional_output", help="If set, will output additional information.", action="store_true", default=False)
    parser.add_argument("--single_snp_model", action="store_true", help="Models are comprised of a single snp per gene", default=False)
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)

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

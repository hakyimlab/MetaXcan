#! /usr/bin/env python
import os
import logging
from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities
from metax.cross_model import Utilities as CrossModelUtilities
from metax.cross_model import JointAnalysis

def run(args):
    if os.path.exists(args.output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.output)
        return
    logging.info("Creating context")
    context = CrossModelUtilities.context_from_args(args)
    genes = context.get_genes()
    results = []

    logging.info("Processing")
    for gene in genes:
        result = JointAnalysis.joint_analysis(context, gene)
        results.append(result)

    results = JointAnalysis.format_results(results)
    Utilities.ensure_requisite_folders(args.output)
    results.to_csv(args.output, index=False, sep="\t")
    logging.info("Ran.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='CrossModel.py %s:'
        'Assess joint analysis of MetaXcan Results' % (__version__))

    parser.add_argument("--cutoff_threshold", help="threshold of variance when truncating SVD (", default=None)
    parser.add_argument("--cutoff_ratio", help="ratio to use when truncating SVD", default=None)
    parser.add_argument("--metaxcan_folder", help="path to metaxcan files", default=None)
    parser.add_argument("--metaxcan_filter", help="regular expression to filter results files", default=[".*csv"], type=str, nargs='+')
    parser.add_argument("--model_product", help="path to file with model feature product", default=None)
    parser.add_argument("--output", help="where you want the output", default=None)
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)

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
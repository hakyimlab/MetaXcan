#! /usr/bin/env python
import os
import logging
from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities
from metax.misc import FeatureMatrix

def _run(args, subset=None, append=None):
    logging.info("Loading expressions")
    manager = FeatureMatrix.build_manager(args.expression_folder, filters = args.expression_filters, standardize=True, subset=subset)

    logging.info("Saving")
    Utilities.ensure_requisite_folders(args.output)
    manager.save_covariances(args.output, append=append)

    logging.info("Ran.")

def run(args):
    if os.path.exists(args.output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.output)
        return

    if args.column_chunks:
        features = sorted(FeatureMatrix.features_in_folder(args.expression_folder))
        total = len(features)
        chunk_size = len(features)/args.column_chunks
        features = [features[i : i+chunk_size] for i in range(0, len(features), chunk_size)]
        l=0
        for i,f in enumerate(features):
            l+=len(f)
            logging.info("Pass %i, %i/%i",i,l,total)
            _run(args, subset=f, append=(i>0))
    else:
        _run(args)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='BuildExpressionProducts.py %s:'
        'Will take "gene expression matrix" files, and build their products (or covariances) with each file taken as a feature' % (__version__))

    parser.add_argument("--column_chunks", help="wether to figure it out on a subset of columsn a t a time (memory reasons)", default=None, type=int)
    parser.add_argument("--expression_folder", help="path to folder with expression", default=None)
    parser.add_argument("--expression_filters", type=str, nargs="+", help="patterns to match expression files", default=["TW_*"])
    parser.add_argument("--output",  help="where you want the output")
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
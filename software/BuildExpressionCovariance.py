#! /usr/bin/env python
import os
import logging
import pandas
import gzip
from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities
from metax.misc import FeatureMatrix

def run(args):
    logging.info("Loading expressions")
    manager = FeatureMatrix.build_manager(args.expression_folder)

    logging.info("Saving")
    manager.save(args.output)

    logging.info("Ran.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='BuildExpressionProducts.py %s:'
        'Will take "gene expression matrix" files, and build their products with each file being a feature' % (__version__))

    parser.add_argument("--expression_folder", help="path to folder with expression", default=None)
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
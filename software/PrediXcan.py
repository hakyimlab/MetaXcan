import logging
from timeit import default_timer as timer

import metax
from metax import Logging
from metax import Exceptions

import Predict
import PrediXcanAssociation

def run(args):
    start = timer()
    logging.info("Running prediction")
    prediction_results = Predict.run(args)

    logging.info("Running association")
    PrediXcanAssociation.run(args, prediction_results)

    end = timer()
    logging.info("Ran PrediXcan in %s seconds"%(str(end-start)))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='PrediXcan.py %s: Single-Tissue PrediXcan' % (metax.__version__))

    Predict.add_arguments(parser)
    PrediXcanAssociation.add_arguments(parser)

    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default=10)
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
#! /usr/bin/env python
import os
import logging
from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import PredictionModel




def run(args):
    if os.path.exists(args.output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.output)
        return

    logging.info("Loading models...")
    model_manager = PredictionModel.load_model_manager(args.models_folder)
    print len(model_manager.get_genes())
    print len(model_manager.get_snps())


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='CovarianceBuilder.py %s:'
        'Collect and process covariance of genotypes' % (__version__))

    parser.add_argument("--models_folder", help="Path to folder with prediction models")
    parser.add_argument("--gtex_genotype_file", help="Path to gtex genotype file")
    parser.add_argument("--gtex_snp_file", help="Path to snp annotation file")
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
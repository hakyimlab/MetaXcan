#! /usr/bin/env python
import os
import logging

from timeit import default_timer as timer

from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities

from metax.predixcan import MultiPrediXcanAssociation
from metax.predixcan import Utilities as MultiPrediXcanUtilities

def run(args):
    start = timer()
    if os.path.exists(args.output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.output)
        return

    with MultiPrediXcanUtilities.mp_context_from_args(args) as context:
        genes = context.get_genes()
        n_genes = len(genes)
        reporter = Utilities.PercentReporter(logging.INFO, n_genes)
        reporter.update(0, "%d %% of model's genes processed so far", force=True)
        results = []
        for i,gene in enumerate(genes):
            r = MultiPrediXcanAssociation.multi_predixcan_association(gene, context)
            results.append(r)
            reporter.update(i, "%d %% of model's genes processed so far")

        results = MultiPrediXcanAssociation.dataframe_from_results(results)

        results.to_csv(args.output, index=False, sep="\t")

    end = timer()
    logging.info("Ran multi tissue predixcan in %s seconds" % (str(end - start)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MTPrediXcan.py %s:'
        'Multi Tissue PrediXcan' % (__version__))

    parser.add_argument("--hdf5_expression_folder", help="Folder with predicted gene expressions.")
    parser.add_argument("--input_phenos_file", help="A text file where on column will be used as phenotype")
    parser.add_argument("--input_phenos_column", help="Name of column from input file to be used as phenotype")
    parser.add_argument("--output", help="File where stuff will be saved.")
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--mode", help="Type of regression. Can be: {}".format(MultiPrediXcanAssociation.MTPMode.K_MODES), default=MultiPrediXcanAssociation.MTPMode.K_LINEAR)

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
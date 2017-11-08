#! /usr/bin/env python
import os
import logging

from timeit import default_timer as timer

from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities

from metax.predixcan import PrediXcanAssociation
from metax.predixcan import Utilities as PrediXcanUtilities

def run(args):
    start = timer()
    if os.path.exists(args.output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.output)
        return

    if (args.hdf5_expression_file and args.expression_file) or \
        (not args.hdf5_expression_file and not args.expression_file):
        logging.info("Provide either hdf5 expression file or plain text expression file")
        return

    with PrediXcanUtilities.p_context_from_args(args) as context:
        genes = context.get_genes()
        n_genes = len(genes)
        reporter = Utilities.PercentReporter(logging.INFO, n_genes)
        reporter.update(0, "%d %% of model's genes processed so far", force=True)
        results = []
        for i,gene in enumerate(genes):
            logging.log(7, "Processing gene %s", gene)
            r = PrediXcanAssociation.predixcan_association(gene, context)
            results.append(r)
            reporter.update(i, "%d %% of model's genes processed so far")
        reporter.update(i, "%d %% of model's genes processed so far")
        results = PrediXcanAssociation.dataframe_from_results(results)
        results = results.fillna("NA")
        results = results.sort_values(by="pvalue")

        Utilities.ensure_requisite_folders(args.output)
        results.to_csv(args.output, index=False, sep="\t", quotechar='"')

    end = timer()
    logging.info("Ran multi tissue predixcan in %s seconds" % (str(end - start)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MTPrediXcan.py %s:'
        'Multi Tissue PrediXcan' % (__version__))

    parser.add_argument("--hdf5_expression_file", help="Folder with predicted gene expressions.")
    parser.add_argument("--expression_file", help="Folder with predicted gene expressions.")
    parser.add_argument("--input_phenos_file", help="A text file where on column will be used as phenotype")
    parser.add_argument('--covariates', type=str, nargs='+',help='Names of covariates in the file', default=[])
    parser.add_argument('--covariates_file', help="File with covariate data. If provided, will force OLS regression")
    parser.add_argument("--input_phenos_column", help="Name of column from input file to be used as phenotype")
    parser.add_argument("--output", help="File where stuff will be saved.")
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--mode", help="Type of regression. Can be: {}".format(PrediXcanAssociation.PMode.K_MODES), default=PrediXcanAssociation.PMode.K_LINEAR)

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
#! /usr/bin/env python
import os
import logging
import traceback

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

    if (args.hdf5_expression_folder and args.expression_folder) or \
        (not args.hdf5_expression_folder and not args.expression_folder):
        logging.info("Provide either hdf5 expression folder or plain text expression folder")
        return

    with MultiPrediXcanUtilities.mp_context_from_args(args) as context:
        genes = context.get_genes()
        n_genes = len(genes)
        reporter = Utilities.PercentReporter(logging.INFO, n_genes)
        reporter.update(0, "%d %% of model's genes processed so far", force=True)

        results = []
        callbacks = {}
        if args.coefficient_output:
            callbacks["coefficient"] = MultiPrediXcanAssociation.SaveCoefs()
        if args.loadings_output:
            callbacks["loadings"] = MultiPrediXcanAssociation.SaveLoadings()

        for i,gene in enumerate(genes):
            if args.MAX_M and i>args.MAX_M-1:
                logging.info("Early exit")
                break

            logging.log(7, "Processing gene %i/%i: %s", i + 1, n_genes, gene)
            if args.throw:
                r = MultiPrediXcanAssociation.multi_predixcan_association(gene, context, callbacks.values())
                results.append(r)
            else:
                try:
                    r = MultiPrediXcanAssociation.multi_predixcan_association(gene, context, callbacks.values())
                    results.append(r)
                except Exception as e:
                    logging.info("Error in gene %s\n%s", gene, traceback.format_exc())
            reporter.update(i, "%d %% of model's genes processed so far")


        reporter.update(i, "%d %% of model's genes processed so far")
        results = MultiPrediXcanAssociation.dataframe_from_results(results, context)
        results = results.sort_values(by="pvalue")
        results = results.fillna("NA")

        Utilities.save_dataframe(results, args.output)
        if args.coefficient_output:
            Utilities.save_dataframe(callbacks["coefficient"].get(), args.coefficient_output)
        if args.loadings_output:
            Utilities.save_dataframe(callbacks["loadings"].get(), args.loadings_output)

    end = timer()
    logging.info("Ran multi tissue predixcan in %s seconds" % (str(end - start)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MTPrediXcan.py %s:'
        'Multi Tissue PrediXcan' % (__version__))

    parser.add_argument("--hdf5_expression_folder", help="Folder with predicted gene expressions. (HDF5 format)")
    parser.add_argument("--expression_folder", help="Folder with predicted gene expressions. (plain text file format)")
    parser.add_argument("--memory_efficient", help="If using plain text expression files, be memory efficient about it. Will be slower.", action="store_true", default=False)
    parser.add_argument("--standardize_expression", help="Standardise input predicted expressions.", action="store_true", default=False)
    parser.add_argument("--code_999", help="values of -999 in expression are to be ignored", action="store_true", default=False)
    parser.add_argument("--expression_pattern", help="Patterns to select expression files", default=None)
    parser.add_argument('--covariates', type=str, nargs='+',help='Names of covariates in the file', default=[])
    parser.add_argument('--covariates_file', help="Text file (or gzip-compressed) with covariate data. If provided, will force OLS regression")
    parser.add_argument("--input_phenos_file", help="Text file (or gzip-compressed) where one column will be used as phenotype")
    parser.add_argument("--input_phenos_column", help="Name of column from input file to be used as phenotype")
    parser.add_argument("--input_phenos_na_values", help="scalar value to be interpreted as 'NA' in the phenotype", nargs="+")
    parser.add_argument("--output", help="File where stuff will be saved.")
    parser.add_argument("--coefficient_output", help="File where coefficients will be saved (when using regularization, built from loadings)")
    parser.add_argument("--loadings_output",help="File where PC loadings will be saved")
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--mode", help="Type of regression. Can be: {}".format(MultiPrediXcanAssociation.MTPMode.K_MODES), default=MultiPrediXcanAssociation.MTPMode.K_LINEAR)
    parser.add_argument("--pc_condition_number", help="Principal components condition number", type=int)
    parser.add_argument("--pc_eigen_ratio", help="Principal components filter, cutoff at proportion to max eigenvalue", type=float)
    parser.add_argument("--MAX_M", help="Compute only the first M genes", type=int)

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
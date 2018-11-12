#! /usr/bin/env python
import os
import logging
import numpy
from timeit import default_timer as timer

import pandas

from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities

from metax.predixcan import MultiPrediXcanAssociation
from metax.predixcan import Utilities as MultiPrediXcanUtilities, Simulations as MultiPredixcanSimulations

def run(args):
    start = timer()

    folder, prefix = os.path.split(args.output_prefix)
    results_name = args.output_prefix + "__mt_results.txt"
    predixcan_results_name = args.output_prefix + "__p_results.txt"
    additional_name = args.output_prefix + "__additional.txt"

    if os.path.exists(results_name):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", results_name)
        return

    #for reproducibility
    numpy.random.seed(100)

    results = []
    additional = []
    predixcan_results = []

    n_max = args.max_n_results
    logging.info("Acquiring context")
    with MultiPredixcanSimulations.context_from_args(args) as context:
        logging.info("processing")
        _c, _cp, _e = context.get_mp_simulation(None)
        for i, gene in enumerate(context.get_genes()):
            if n_max and i+1>n_max:
                logging.info("Max runs met")
                break
            logging.log(9, "%d Gene %s", i, gene)
            r, add, p = MultiPredixcanSimulations.simulate(gene, context)
            if r is None:
                logging.log(9, "%s could not be simulated", gene)
                continue
            results.append(r)
            additional.append(add)

            if p is not None:
                predixcan_results.append(p)

    results = MultiPrediXcanAssociation.dataframe_from_results(results, _c).sort_values(by="pvalue")
    additional = pandas.concat(additional)

    Utilities.ensure_requisite_folders(results_name)
    Utilities.save_dataframe(results, results_name)
    Utilities.save_dataframe(additional, additional_name)

    if len(predixcan_results):
        predixcan_results = pandas.concat(predixcan_results)
        Utilities.save_dataframe(predixcan_results, predixcan_results_name)
    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='MPSimulation.py')

    parser.add_argument("--expression_folder", help="Folder with predicted gene expressions. (plain text file format)")
    parser.add_argument("--expression_pattern", help="Patterns to select expression files", default=None)
    parser.add_argument("--input_phenos_file", help="Text file (or gzip-compressed) where one column will be used as phenotype")
    parser.add_argument("--simulation_type", help="What kind of genotype to simulate: [random, combination, simple]")
    parser.add_argument("--output_prefix", help="File where stuff will be saved.")
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--code_999", help="values of -999 in expression are to be ignored", action="store_true", default=False)
    parser.add_argument("--mode", help="Type of regression. Can be: {}".format(MultiPrediXcanAssociation.MTPMode.K_MODES), default=MultiPrediXcanAssociation.MTPMode.K_LINEAR)
    parser.add_argument("--pc_condition_number", help="Principal components condition number", type=int)
    parser.add_argument("--pc_eigen_ratio", help="Principal components filter, cutoff at proportion to max eigenvalue", type=float)
    parser.add_argument("--standardize_expression", help="Standardise input predicted expressions.", action="store_true", default=False)
    parser.add_argument("--only_truth", help="Run Multi-PrediXcan only with selected causal models.", action="store_true", default=False)
    parser.add_argument("--simulation_parameters", help="Depends on particular scheme", action="append", nargs=2)
    parser.add_argument("--do_predixcan", help="Also compute predixcan association", action="store_true", default=False)
    parser.add_argument("--max_n_results", help="Optional. If provided, run up to as many analysis", type=int)

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
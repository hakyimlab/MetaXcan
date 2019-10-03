#! /usr/bin/env python
import os
import logging
import traceback

from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import Utilities
from metax.gwas import Utilities as GWASUtilities
from metax.cross_model import Utilities as CrossModelUtilities
from metax.cross_model import JointAnalysis

from timeit import default_timer as timer

def run(args):
    start = timer()
    if os.path.exists(args.output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.output)
        return
    logging.info("Creating context")
    context = CrossModelUtilities.context_from_args(args)
    results = []

    n_genes = context.get_n_genes()
    reporter = Utilities.PercentReporter(logging.INFO, n_genes)

    logging.info("Processing")

    reporter.update(0, "%d %% of model's genes processed so far")
    for i,gene in enumerate(context.get_genes()):
        if args.MAX_M and i > args.MAX_M-1:
            logging.info("Early exit")
            break

        logging.log(7, "Gene %d/%d: %s", i+1, n_genes, gene)
        if args.throw:
            result = JointAnalysis.joint_analysis(context, gene)
            results.append(result)
        else:
            try:
                result = JointAnalysis.joint_analysis(context, gene)
                results.append(result)
            except Exception as e:
                logging.info("Error in gene %s\n%s", gene, traceback.format_exc())
        reporter.update(i, "%d %% of model's genes processed so far")

    results = JointAnalysis.format_results(results)
    Utilities.ensure_requisite_folders(args.output)
    results.to_csv(args.output, index=False, sep="\t", na_rep="NA")

    end = timer()
    logging.info("Ran multi tissue in %s seconds" % (str(end - start)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='CrossModel.py %s:'
        'Assess joint analysis of MetaXcan Results' % (__version__))

    parser.add_argument("--cleared_snps", help="SNPS to analyze. If you don't provide this, you must provide GWAS and models.")
    parser.add_argument("--models_folder", help="Path to folder with prediction models")
    parser.add_argument("--models_name_filter", help="Path to folder with prediction models", type=str, nargs='+')
    parser.add_argument("--models_name_pattern", help="regular expression to detect tissue name from file names")
    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")

    parser.add_argument("--gwas_folder", help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.")
    parser.add_argument("--gwas_file_pattern", help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).")
    parser.add_argument("--gwas_file", help="Path to GWAS file; alternative to gwas_folder and gwas_pattern")
    GWASUtilities.add_gwas_arguments_to_parser(parser)

    parser.add_argument("--regularization", help="Add a regularization term to correct for expression covariance matrix singularity", default=None, type=float)
    parser.add_argument("--cutoff_condition_number", help="condition number of eigen values to use when truncating SVD", default=None, type=float)
    parser.add_argument("--cutoff_eigen_ratio", help="ratio to use when truncating SVD", default=None, type=float)
    parser.add_argument("--cutoff_threshold", help="threshold of variance when truncating SVD", default=None, type=float)
    parser.add_argument("--cutoff_trace_ratio", help="ratio to use when truncating SVD", default=None, type=float)
    parser.add_argument("--metaxcan_folder", help="path to metaxcan files", default=None)
    parser.add_argument("--metaxcan_filter", help="regular expression to filter results files", default=[".*csv"], type=str, nargs='+')
    parser.add_argument("--metaxcan_file_name_parse_pattern", help="Optional regular expression to get phenotype name and model name from MetaXcan result files.", default = None)
    parser.add_argument("--model_product", help="path to file with model feature product", default=None)
    parser.add_argument("--permissive_model_product", action="store_true", help="Tells the Model Product to need some slack, entries may be missing on numerical problems", default=False)
    parser.add_argument("--snp_covariance", help="path  to snp covariance", default=None)
    parser.add_argument("--output", help="where you want the output", default=None)
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--throw", action="store_true", help="Throw exception on error", default=False)
    parser.add_argument("--trimmed_ensemble_id", action="store_true", help="Use ensemble ids without version", default=False)
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
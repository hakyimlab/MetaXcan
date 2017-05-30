#! /usr/bin/env python
import os
import logging

from timeit import default_timer as timer

from metax import __version__
from metax import Logging
from metax import Exceptions
from metax import PredictionModel
from metax import Utilities
from metax import MatrixManager
from metax.genotype import  Utilities as GenotypeUtilities
from metax.genotype import GenotypeAnalysis

def run(args):
    if os.path.exists(args.snp_covariance_output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.snp_covariance_output)
        return

    if os.path.exists(args.gene_variance_output):
        logging.info("%s already exists, you have to move it or delete it if you want it done again", args.gene_variance_output)
        return

    start = timer()

    logging.info("Loading models...")
    model_manager = PredictionModel.load_model_manager(args.models_folder)
    all_snps = model_manager.get_rsids()

    logging.info("processing genotype")
    for chromosome, metadata, dosage in GenotypeUtilities.genotype_by_chromosome_from_args(args, all_snps):
        logging.log(9, "Processing chromosome %s", str(chromosome))
        covariance_results = []
        variance_results = []

        context = GenotypeAnalysis.GenotypeAnalysisContext(metadata, dosage, model_manager)
        genes = context.get_genes()
        reporter = Utilities.PercentReporter(9, len(genes))
        reporter.update(0, "%d %% of genes processed so far in chromosome " + str(chromosome))
        for i,gene in enumerate(genes):
            reporter.update(i, "%d %% of genes processed so far in chromosome "+str(chromosome))
            cov_data = GenotypeAnalysis.get_prediction_covariance(context, gene)
            cov_data = MatrixManager._flatten_matrix_data([cov_data])
            covariance_results.extend(cov_data)

            variance_data = GenotypeAnalysis.get_prediction_variance(context, gene)
            variance_results.extend(variance_data)
        reporter.update(len(genes), "%d %% of genes processed so far in chromosome " + str(chromosome))

        logging.log(9, "writing chromosome results")

        Utilities.save_table(covariance_results, args.snp_covariance_output,
                                    mode="w" if chromosome ==1 else "a",
                                    header= GenotypeAnalysis.COVARIANCE_COLUMNS if chromosome==1 else None)

        Utilities.save_table(variance_results, args.gene_variance_output,
                                    mode="w" if chromosome == 1 else "a",
                                    header=GenotypeAnalysis.VARIANCE_COLUMNS if chromosome == 1 else None)

    end = timer()
    logging.info("Ran covariance builder in %s seconds" % (str(end - start)))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='CovarianceBuilder.py %s:'
        'Collect and process covariance of genotypes' % (__version__))

    parser.add_argument("--models_folder", help="Path to folder with prediction models")
    parser.add_argument("--gtex_genotype_file", help="Path to gtex genotype file")
    parser.add_argument("--gtex_snp_file", help="Path to snp annotation file")
    parser.add_argument("--snp_covariance_output", help="where you want the output", default=None)
    parser.add_argument("--gene_variance_output", help="where you want the output", default=None)
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
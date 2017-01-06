#! /usr/bin/env python
__author__ = 'heroico'
import metax
__version__ = metax.__version__

import logging
import pandas
import os

from metax import PredictionModel
from metax import MatrixManager
from metax import Logging
from metax import Utilities
from metax import Exceptions
from metax.metaxcan import AssociationCalculation

class MResult(object):
    def __init__(self):
        self.gene = None,
        self.gene_name = None
        self.zscore = None
        self.effect_size = None
        self.p = None
        self.VAR_g = None
        self.n = None
        self.n_cov = None
        self.n_model = None
        self.gene_R2 = None
        self.gene_p = None
        self.gene_q = None

    HEADER="gene,gene_name,zscore,effect_size,pvalue,VAR_g,pred_perf_R2,pred_perf_p,pred_perf_q,n_snps_used,n_snps_in_cov,n_snps_in_model\n"

    def toCSVLine(self):
        line = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
               (self.gene, self.gene_name, self.zscore, self.effect_size, self.p, self.VAR_g, self.gene_R2, self.gene_p, self.gene_q, self.n, self.n_cov, self.n_model)
        return line

def _beta_loader(args):
    beta_contents = Utilities.contentsWithPatternsFromFolder(args.beta_folder, [])
    for beta_name in beta_contents:
        logging.info("Processing %s", beta_name)
        beta_path = os.path.join(args.beta_folder, beta_name)
        b = pandas.read_table(beta_path)
        yield b

def _gwas_wrapper(gwas):
    logging.info("Processing input gwas")
    yield gwas

def run(args, _gwas=None):
    if os.path.exists(args.output_file):
        logging.info("%s already exists, move it or delete it if you want it done again")
        return
    logging.info("Started metaxcan association")

    logging.info("Loading model from: %s", args.model_db_path)
    model = PredictionModel.load_model(args.model_db_path)

    logging.info("Loading covariance data from: %s", args.covariance)
    covariance_manager = MatrixManager.MatrixManager(args.covariance)

    total_snps = len(set(model.weights.rsid))
    snps_found=set()
    reporter = Utilities.PercentReporter(logging.INFO, total_snps)

    results = []
    gwas_gen = _gwas_wrapper(_gwas) if _gwas is not None else _beta_loader(args)
    for gwas in gwas_gen:
        context = AssociationCalculation.Context(gwas, model, covariance_manager)
        i_genes, i_snps = AssociationCalculation.intersection_d(model, gwas)
        snps_found.update(i_snps)
        for gene in i_genes:
            r = AssociationCalculation.association(gene, context)
            results.append(r)
        reporter.update(len(snps_found), "%d %% of model's snps found so far in the gwas study")
    results = AssociationCalculation.dataframe_from_results(results)
    results.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M04_zscores.py %s: Build ZScores from GWAS data.' % (__version__,))

    parser.add_argument("--model_db_path",
                        help="name of weight db in data folder",
                        default=None)

    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default=None)

    parser.add_argument("--beta_folder",
                        help="name of folder containing beta data",
                        default=None)

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")

    parser.add_argument("--zscore_scheme",
                        help="Scheme for zscore calculation. Options are:"
                             "'beta_z' (uses zscore of beta and sigma_l from input file), default;"
                            " 'beta_z_and_ref' (uses zscore of beta and sigma_l from reference population);"
                            " 'metaxcan' ('bare metal' MetaXcan, normalization recommended); "
                            " 'metaxcan_from_reference' ('bare metal' MetaXcan, using reference variance);",
                        default=None)

    parser.add_argument("--normalization_scheme",
                        help="Scheme for zscore normalization, relevant for 'global normalization' scheme. Options are:"
                            "'none';"
                            " 'from_pheno', estimate normalization constant from phenotype file, needs 'sigma_l' and 'standard error' in phenotype;"
                            " 'from_reference', estimate normalization constant from reference, needs 'standard error' on phenotype",
                        default=None)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--keep_ens_version",
                        help="If set, will keep the -version- postfix in gene id.",
                    action="store_true",
                    default=False)

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)


    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))
    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException as e:
            logging.error("Error:%s", e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)

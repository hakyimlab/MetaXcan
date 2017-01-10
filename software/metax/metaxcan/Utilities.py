import logging
import pandas
import os

import AssociationCalculation

from .. import Constants
from .. import Utilities
from .. import MatrixManager
from ..PredictionModel import WDBQF, load_model

class SimpleContext(object):
    def __init__(self, gwas, model, covariance):
        self.gwas = gwas
        self.model = model
        self.covariance = covariance

    def get_weights(self, gene):
        w = self.model.weights
        w = w[w.gene == gene]
        return w

    def get_covariance(self, gene, snps):
        return self.covariance.get(gene, snps)

    def get_n_in_covariance(self, gene):
        return self.covariance.n_snps(gene)

    def get_gwas(self, snps):
        g = self.gwas
        g = g[g[Constants.SNP].isin(snps)]
        return g

    def get_model_snps(self):
        return set(self.model.weights.rsid)

    def get_data_intersection(self):
        weights, gwas = self.model.weights, self.gwas
        k = pandas.merge(weights, gwas, how='inner', left_on="rsid", right_on="snp")
        genes = k.gene.drop_duplicates().values
        snps = k.rsid.drop_duplicates().values
        return genes, snps

def _prepare_gwas(gwas):
    #If zscore is numeric, then everything is fine with us.
    # if not, try to remove "NA" strings.
    try:
        gwas = gwas[gwas.zscore != "NA"]
    except Exception as e:
        logging.log(9, "Unexpected issue preparing gwas... %s", str(e))
        pass
    return gwas

def _prepare_model(model):
    K = WDBQF.K_GENE
    g = model.weights[K]
    model.weights[K] = pandas.Categorical(g, g.drop_duplicates())
    return model

def _beta_loader(args):
    beta_contents = Utilities.contentsWithPatternsFromFolder(args.beta_folder, [])
    r = pandas.DataFrame()
    for beta_name in beta_contents:
        logging.info("Processing %s", beta_name)
        beta_path = os.path.join(args.beta_folder, beta_name)
        b = pandas.read_table(beta_path)
        r = pandas.concat([r, b])
    return r

def _gwas_wrapper(gwas):
    logging.info("Processing input gwas")
    return gwas

def build_context(args, gwas):
    logging.info("Loading model from: %s", args.model_db_path)
    model = load_model(args.model_db_path)

    logging.info("Loading covariance data from: %s", args.covariance)
    covariance_manager = MatrixManager.load_matrix_manager(args.covariance)

    gwas = _gwas_wrapper(gwas) if gwas is not None else _beta_loader(args)
    context = _build_context(model, covariance_manager, gwas)
    return context

def _build_context(model, covariance_manager, gwas):
    model = _prepare_model(model)
    gwas = _prepare_gwas(gwas)
    context = SimpleContext(gwas, model, covariance_manager)
    return context



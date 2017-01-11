import logging
import pandas
import os
import numpy

from .. import Constants
from .. import Utilities
from .. import MatrixManager
from ..PredictionModel import WDBQF, load_model, dataframe_from_weight_data

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
        return self.covariance.get(gene, snps, strict=False)

    def get_n_in_covariance(self, gene):
        return self.covariance.n_snps(gene)

    def get_gwas(self, snps):
        g = self.gwas
        g = g[g[Constants.SNP].isin(snps)]
        return g

    def get_model_snps(self):
        return set(self.model.weights.rsid)

    def get_data_intersection(self):
        return _data_intersection(self.model, self.gwas)

    def provide_calculation(self, gene):
        w = self.get_weights(gene)
        gwas = self.get_gwas(w[WDBQF.K_RSID].values)

        i = pandas.merge(w, gwas, left_on="rsid", right_on="snp")
        if not Constants.BETA in i: i[Constants.BETA] = None
        i = i[[Constants.SNP, WDBQF.K_WEIGHT, Constants.ZSCORE, Constants.BETA]]

        snps, cov = self.get_covariance(gene, i[Constants.SNP].values)

        # fast subsetting and aligning
        d_columns = i.columns.values
        if snps is not None and len(snps):
            d = {x[0]: x for x in i.values}
            d = [d[snp] for snp in snps]
            d = zip(*d)
            d = {d_columns[i]:d[i] for i in xrange(0, len(d_columns))}
            i = pandas.DataFrame(d)
        else:
            i = pandas.DataFrame(columns=d_columns)
        return w, i, cov, snps

class OptimizedContext(SimpleContext):
    def __init__(self, gwas, model, covariance):
        self.covariance = covariance
        self.weight_data, self.snps_in_model = _prepare_weight_data(model)
        self.gwas_data = _prepare_gwas_data(gwas)

    def get_weights(self, gene):
        w = self.weight_data[gene]
        w = dataframe_from_weight_data(zip(*w))
        return w

    def get_model_snps(self):
        return set(self.snps_in_model)

    def get_gwas(self, snps):
        snps = set(snps)
        g = self.gwas_data
        g = [g[x] for x in snps if x in g]
        if len(g):
            g = zip(*g)
            g = pandas.DataFrame({Constants.SNP:g[0], Constants.ZSCORE:g[1], Constants.BETA:g[2]})
        else:
            g = pandas.DataFrame(columns=[Constants.SNP, Constants.ZSCORE, Constants.BETA])
        return g

    def get_data_intersection(self):
        return _data_intersection_2(self.weight_data, self.gwas_data)

def _data_intersection(model, gwas):
    weights = model.weights
    k = pandas.merge(weights, gwas, how='inner', left_on="rsid", right_on="snp")
    genes = k.gene.drop_duplicates().values
    snps = k.rsid.drop_duplicates().values
    return genes, snps

def _data_intersection_2(weight_data, gwas_data):
    genes = set()
    snps = set()
    for gene, entries in weight_data.iteritems():
        gs = zip(*entries)[WDBQF.RSID]
        for s in gs:
            if s in gwas_data:
                genes.add(gene)
                snps.add(s)
    return genes, snps

def _prepare_gwas(gwas):
    #If zscore is numeric, then everything is fine with us.
    # if not, try to remove "NA" strings.
    try:
        gwas = gwas[gwas.zscore != "NA"]
    except Exception as e:
        logging.log(9, "Unexpected issue preparing gwas... %s", str(e))
        pass

    if not Constants.BETA in gwas:
        gwas[Constants.BETA] = numpy.nan

    return gwas

def _prepare_gwas_data(gwas):
    gwas = gwas[[Constants.SNP, Constants.ZSCORE, Constants.BETA]]
    data = {}
    for x in gwas.values:
        data[x[0]] = x
    return data

def _prepare_model(model):
    K = WDBQF.K_GENE
    g = model.weights[K]
    model.weights[K] = pandas.Categorical(g, g.drop_duplicates())
    return model

def _prepare_weight_data(model):
    d = {}
    snps = set()
    for x in model.weights.values:
        gene = x[WDBQF.GENE]
        if not gene in d:
            d[gene] = []
        entries = d[gene]
        entries.append(x)
        snps.add(x[WDBQF.RSID])
    return d, snps

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
    gwas = _prepare_gwas(gwas)
    context = OptimizedContext(gwas, model, covariance_manager)
    return context

def _build_simple_context(model, covariance_manager, gwas):
    model = _prepare_model(model)
    gwas = _prepare_gwas(gwas)
    context = SimpleContext(gwas, model, covariance_manager)
    return context


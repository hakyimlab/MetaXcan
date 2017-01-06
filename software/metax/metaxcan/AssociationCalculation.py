import logging
import pandas
import numpy
from numpy import dot as d

from ..PredictionModel import WDBEQF, WDBQF
from .. import Constants

class ARF(object):
    """Association result format"""
    GENE = 0
    ZSCORE = 1
    EFFECT_SIZE = 2
    N_SNPS_IN_MODEL = 3
    N_SNPS_IN_COV = 4
    N_SNPS_USED = 5

    K_GENE = "gene"
    K_ZSCORE = "zscore"
    K_EFFECT_SIZE = "effect_size"
    K_N_SNPS_IN_MODEL = "n_snps_in_model"
    K_N_SNPS_IN_COV = "n_snps_in_cov"
    K_N_SNPS_USED = "n_snps_used"

    order=[(GENE,K_GENE), (ZSCORE,K_ZSCORE), (EFFECT_SIZE,K_EFFECT_SIZE), (N_SNPS_IN_MODEL, K_N_SNPS_IN_MODEL), (N_SNPS_IN_COV, K_N_SNPS_IN_COV), (N_SNPS_USED,K_N_SNPS_USED)]


class Context(object):
    def __init__(self, gwas=None, prediction_model=None, covariance_manager=None):
        self.gwas = gwas
        self.prediction_model = prediction_model
        self.covariance_manager = covariance_manager

def intersection_d(prediction_model, gwas):
    k = pandas.merge(prediction_model.weights, gwas, how='inner', left_on="rsid", right_on="snp")
    genes = k.gene.drop_duplicates().values
    snps = k.rsid.drop_duplicates().values
    return genes, snps

def association(gene, context):
    #capture context
    model = context.prediction_model
    covariance = context.covariance_manager
    gwas = context.gwas

    # Select and align data for gene
    w = model.weights[model.weights.gene == gene]
    i = pandas.merge(w, gwas, left_on="rsid", right_on="snp")
    if not Constants.BETA in i: i[Constants.BETA] = None
    i = i[[Constants.SNP, WDBQF.K_WEIGHT, Constants.ZSCORE, Constants.BETA]]
    i = i.dropna()
    i = i[i.zscore != "NA"]
    snps, cov = covariance.get(gene, i[Constants.SNP].values)
    i = i[i.snp.isin(snps)]
    i.snp = pandas.Categorical(i.snp, snps)
    i = i.sort_values(by='snp')

    #some stats
    n_snps_in_model = len(w.effect_allele)
    n_snps_used = len(i.weight)
    n_snps_in_cov = covariance.n_snps(gene)

    zscore = None
    effect_size = None
    sigma_g_2 = None
    if n_snps_used > 0:
        # get sigma from reference
        variances = numpy.diag(cov)
        i['sigma_l'] = numpy.sqrt(variances)

        #da calculeishon
        sigma_g_2 = float(d(d(i.weight,cov),i.weight))
        if sigma_g_2 >0:
            try:
                zscore = numpy.sum(i.weight * i.zscore * i.sigma_l) / numpy.sqrt(sigma_g_2)
                effect_size = numpy.sum(i.weight * i.beta * (i.sigma_l**2))/ sigma_g_2
            except Exception as e:
                logging.log(9, "Unexpected exception when calculating zscore: %s, %s", gene, str(e))

    return gene, zscore, effect_size, n_snps_in_model, n_snps_in_cov, n_snps_used

def dataframe_from_results(results):
    r = pandas.DataFrame({key: results[order] for order, key in ARF.order})
    r = r[[key for order,key in ARF.order]]
    return r
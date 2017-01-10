import logging
import pandas
import numpy
from numpy import dot as d

from .. import Constants
from .. import Exceptions

from ..PredictionModel import WDBEQF, WDBQF

class ARF(object):
    """Association result format"""
    GENE = 0
    ZSCORE = 1
    EFFECT_SIZE = 2
    SIGMA_G_2 = 3
    N_SNPS_IN_MODEL = 4
    N_SNPS_IN_COV = 5
    N_SNPS_USED = 6

    K_GENE = "gene"
    K_ZSCORE = "zscore"
    K_EFFECT_SIZE = "effect_size"
    K_VAR_G = "var_g"
    K_N_SNPS_IN_MODEL = "n_snps_in_model"
    K_N_SNPS_IN_COV = "n_snps_in_cov"
    K_N_SNPS_USED = "n_snps_used"

    order=[(GENE,K_GENE), (ZSCORE,K_ZSCORE), (EFFECT_SIZE,K_EFFECT_SIZE), (SIGMA_G_2, K_VAR_G), (N_SNPS_IN_MODEL, K_N_SNPS_IN_MODEL), (N_SNPS_IN_COV, K_N_SNPS_IN_COV), (N_SNPS_USED,K_N_SNPS_USED)]

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract context")
    def get_weights(self, gene): pass
    def get_covariance(self, gene, snps): pass
    def get_n_in_covariance(self, gene): pass
    def get_gwas(self, snps): pass
    def get_model_snps(self): pass
    def get_data_intersection(self): pass

def association(gene, context, return_snps=False):
    #capture context

    # Select and align data for gene
    w = context.get_weights(gene)
    gwas = context.get_gwas(w[WDBQF.K_RSID].values)
    i = pandas.merge(w, gwas, left_on="rsid", right_on="snp")
    if not Constants.BETA in i: i[Constants.BETA] = None
    i = i[[Constants.SNP, WDBQF.K_WEIGHT, Constants.ZSCORE, Constants.BETA]]

    snps, cov = context.get_covariance(gene, i[Constants.SNP].values)
    if snps is None:
        r = (gene, numpy.nan, numpy.nan, len(w.effect_allele), 0, 0)
        if return_snps:
            return r, set(i.snp)
        else:
            return r

    i = i[i.snp.isin(set(snps))]
    i.snp = pandas.Categorical(i.snp, snps)
    i = i.sort_values(by='snp')

    #some stats
    n_snps_in_model = len(w.effect_allele)
    n_snps_used = len(i.weight)
    n_snps_in_cov = context.get_n_in_covariance(gene)

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

    r = (gene, zscore, effect_size, sigma_g_2, n_snps_in_model, n_snps_in_cov, n_snps_used)
    if return_snps:
        return r, set(i.snp)
    else:
        return r

def dataframe_from_results(results):
    r = pandas.DataFrame({key: results[order] for order, key in ARF.order})
    r = r[[key for order,key in ARF.order]]
    return r
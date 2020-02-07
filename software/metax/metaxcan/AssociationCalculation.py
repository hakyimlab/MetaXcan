import logging
import pandas
import numpy
from numpy import dot as d
from scipy import stats

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

class ASF(object):
    """Additional Stats format"""
    GENE = 0
    BEST_GWAS_P = 1
    LARGEST_WEIGHT = 2

    K_GENE = "gene"
    K_BEST_GWAS_P = "best_gwas_p"
    K_LARGEST_WEIGHT = "largest_weight"

    order=[(GENE,K_GENE), (BEST_GWAS_P, K_BEST_GWAS_P), (LARGEST_WEIGHT,K_LARGEST_WEIGHT)]

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract context")
    def get_weights(self, gene): pass
    def get_covariance(self, gene, snps): pass
    def get_n_in_covariance(self, gene): pass
    def get_gwas(self, snps): pass
    def get_model_snps(self): pass
    def get_data_intersection(self): pass
    def provide_calculation(self, gene): pass
    def get_model_info(self): pass

def association(gene, context, return_snps=False):
    #capture context
    n_snps_in_model, i, cov, snps = context.provide_calculation(gene)

    if logging.getLogger().getEffectiveLevel() < 10:
        d_ = numpy.linalg.eig(cov)[0]
        if numpy.sum(numpy.less(d_,1e-6)):
            logging.info("Gene %s has covariance close to singular", gene)

    #some stats
    snps_used = i[Constants.SNP]
    n_snps_used = len(snps_used)
    n_snps_in_cov = context.get_n_in_covariance(gene)

    zscore, effect_size, sigma_g_2 = numpy.nan, numpy.nan, numpy.nan

    if n_snps_used > 0:
        i_weight = i[WDBQF.K_WEIGHT]
        i_zscore = i[Constants.ZSCORE]
        i_beta = i[Constants.BETA]
        # sigma from reference
        variances = numpy.diag(cov)
        i_sigma_l = numpy.sqrt(variances)

        #da calcooleishon
        sigma_g_2 = float(d(d(i_weight,cov),i_weight))

        if sigma_g_2 >0:
            try:
                zscore = numpy.sum(i_weight * i_zscore * i_sigma_l) / numpy.sqrt(sigma_g_2)
                effect_size = numpy.sum(i_weight * i_beta * (i_sigma_l**2))/ sigma_g_2
            except Exception as e:
                logging.log(9, "Unexpected exception when calculating zscore: %s, %s", gene, str(e))

    r = (gene, zscore, effect_size, sigma_g_2, n_snps_in_model, n_snps_in_cov, n_snps_used)

    if return_snps:
        return r, set(snps_used)
    else:
        return r

def dataframe_from_results(results):
    results = list(zip(*results))
    if len(results) == 0:
        return pandas.DataFrame({key:[] for order,key in ARF.order})

    r = pandas.DataFrame({key: results[order] for order, key in ARF.order})
    r = r[[key for order,key in ARF.order]]
    return r

def additional_stats(gene, context):
    #capture context
    n_snps_in_model, i, cov, snps = context.provide_calculation(gene)

    #some stats
    snps_used = i[Constants.SNP]
    n_snps_used = len(snps_used)

    if n_snps_used == 0:
        return gene, numpy.nan, numpy.nan

    i_zscore = i[Constants.ZSCORE]
    best_zscore = numpy.max(numpy.abs(i_zscore))
    p = 2 * stats.norm.cdf(-best_zscore)

    i_weight = i[WDBQF.K_WEIGHT]
    best_weight = numpy.max(numpy.abs(i_weight))

    return gene, p, best_weight

def dataframe_from_aditional_stats(stats):
    stats = list(zip(*stats))
    if len(stats) == 0:
        return pandas.DataFrame({key:[] for order,key in ASF.order})

    r = pandas.DataFrame({key: stats[order] for order, key in ASF.order})
    r = r[[key for order,key in ASF.order]]
    return r

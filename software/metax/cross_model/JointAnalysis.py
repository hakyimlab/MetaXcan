import numpy
from numpy.core import dot, array
from scipy import stats

from .. import Exceptions
from ..misc import Math
from .. import Utilities

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract Joint Analysis context")
    def get_genes(self): pass
    def get_metaxcan_zscores(self, gene): pass
    def get_model_matrix(self, gene, tissues): pass

class CalculationStatus(object):
    OK=0
    NO_DATA=-1
    NO_METAXCAN_RESULTS=-2
    NO_PRODUCT=-3

def joint_analysis(context, gene):
    pvalue, n, n_indep, p_i_best, p_i_worst, status = None, None, None, None, None, CalculationStatus.NO_DATA
    g = gene.split(".")[0]

    zscores, tissue_labels = context.get_metaxcan_zscores(gene)
    if not zscores or len(zscores) == 0:
        status = CalculationStatus.NO_METAXCAN_RESULTS
        return g, pvalue, n, n_indep, p_i_best, p_i_worst, status
    n = len(zscores)

    labels, matrix = context.get_model_matrix(gene, tissue_labels)
    if not labels or len(labels) == 0:
        status = CalculationStatus.NO_PRODUCT
        return g, pvalue, n, n_indep, p_i_best, p_i_worst, status

    trace = numpy.trace(matrix)
    cutoff = 0.1*trace
    zscores = array([z for i,z in enumerate(zscores) if tissue_labels[i] in labels])
    inv, n_indep, eigen = Math.capinv(matrix, cutoff)

    max_z = numpy.max(numpy.abs(zscores))
    p_i_best = 2*stats.norm.cdf(-max_z)
    min_z = numpy.min(numpy.abs(zscores))
    p_i_worst = 2*stats.norm.cdf(-min_z)

    w = float(dot(dot(zscores, inv), zscores))
    pvalue = 1 - stats.chi2.cdf(w, n_indep)
    status = CalculationStatus.OK

    return g, pvalue, n, n_indep, p_i_best, p_i_worst, status

def format_results(results):
    columns = ["gene", "pvalue", "n", "n_indep", "p_i_best", "p_i_worst", "status"]
    results = Utilities.to_dataframe(results, columns)
    results = results.sort_values(by="pvalue")
    results = results.fillna("NA")
    return results
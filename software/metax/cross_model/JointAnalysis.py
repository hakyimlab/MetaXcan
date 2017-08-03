import numpy
import logging
from numpy.core import dot, array
from scipy import stats

from .. import Exceptions
from ..misc import Math
from .. import Utilities

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract Joint Analysis context")
    def get_genes(self): raise  Exceptions.NotImplemented("Context: get_genes")
    def get_n_genes(self): raise  Exceptions.NotImplemented("Context: get_n_genes")
    def get_metaxcan_zscores(self, gene): raise  Exceptions.NotImplemented("Context: get_metaxcan_zscores")
    def get_model_matrix(self, gene, tissues): raise  Exceptions.NotImplemented("Context: get_model_matrix")
    def get_cutoff(self, matrix): raise  Exceptions.NotImplemented("Context: get_cutoff")

class ContextMixin(object):
    def __init__(self):
        self.metaxcan_results_manager = None
        self.matrix_manager = None
        self.cutoff = None

    def get_metaxcan_zscores(self, gene):
        if "." in gene: gene = gene.split(".")[0]
        results = self.metaxcan_results_manager.results_for_gene(gene)
        return results

    def get_model_matrix(self, gene, tissues):
        return self.matrix_manager.get(gene, tissues)

    def get_cutoff(self, matrix):
        return self.cutoff(matrix)

class CalculationStatus(object):
    OK=0
    NO_DATA=-1
    NO_METAXCAN_RESULTS=-2
    NO_PRODUCT=-3
    INSUFFICIENT_NUMERICAL_RESOLUTION = -4
    SINGULAR_COVARIANCE = -5
    INVERSE_ERROR = -6
    COMPLEX_COVARIANCE = -7
    INADEQUATE_INVERSE = -8


# Todo: remove
DEBUG=False
#TODO remove
import pandas
import matplotlib.pyplot as plt
def get_merged(gene, labels, zscores):
    predixcan = pandas.read_table("data/cross_model_p/implementation_prototype/results_t1d/predixcan.txt")
    p = predixcan[predixcan.gene == gene]
    p.gene = p.gene.str.split(".").str.get(0)
    p.model = p.model.str.split("TW_").str.get(1)
    m = pandas.DataFrame(data={"model":labels, "zscore":zscores})
    merged = p.merge(m, on="model")
    return merged

def plot(merged):
    plt.plot(merged.zscore_x, merged.zscore_y, 'ro')
    plt.plot([numpy.min(merged.zscore_x), numpy.max(merged.zscore_x)],
             [numpy.min(merged.zscore_x), numpy.max(merged.zscore_x)])
    plt.show()

def joint_analysis(context, gene):
    g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status \
        = None, None, None, None, None, None, None, None, None, None, None, None, None, CalculationStatus.NO_DATA
    g = gene.split(".")[0]

    if DEBUG:
        if not g in ["ENSG00000171824"]:
        #if not g in ["ENSG00000053371"]: #this is close, previous is not
            return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status
        print "a"
        from IPython import embed; embed()

    zscores, tissue_labels = context.get_metaxcan_zscores(gene)
    if not zscores or len(zscores) == 0:
        status = CalculationStatus.NO_METAXCAN_RESULTS
        return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status
    n = len(zscores)
    z_min = numpy.min(zscores)
    z_max = numpy.max(zscores)
    z_mean = numpy.mean(zscores)
    if (len(zscores)>1):
        z_sd = numpy.std(zscores, ddof=1)

    labels, matrix = context.get_model_matrix(gene, tissue_labels)

    if not labels or len(labels) == 0:
        status = CalculationStatus.NO_PRODUCT
        return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status

    # also, check that the matrix actually makes sense. We are currently returning it just in case but matrices with complex covariance are suspicious.
    e, v = numpy.linalg.eig(matrix)
    if numpy.imag(e).any():
        status = CalculationStatus.COMPLEX_COVARIANCE
        e = numpy.real(e)
        eigen_max, eigen_min = numpy.max(e), numpy.min(e)
        return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status

    eigen_max, eigen_min = numpy.max(e), numpy.min(e)

    cutoff = context.get_cutoff(matrix)
    _d = {tissue_labels[i]:zscores[i] for i in xrange(0, len(tissue_labels))}
    zscores = array([_d[l] for l in labels])
    inv, n_indep, eigen = Math.capinv(matrix, cutoff, context.epsilon)

    max_z = numpy.max(numpy.abs(zscores))
    p_i_best = 2*stats.norm.cdf(-max_z)
    min_z = numpy.min(numpy.abs(zscores))
    p_i_worst = 2*stats.norm.cdf(-min_z)

    #TODO: implement a better heuristic
    try:
        eigen_w, eigen_v = numpy.linalg.eig(inv)
    except:
        #WTCCC 'ENSG00000204560.5'
        logging.log(8, "Problems with inverse for %s, skipping", gene)
        status = CalculationStatus.INVERSE_ERROR
        return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status

    # if numpy.max(eigen_w) > 1e10:
    #     logging.log(8, "gene %s has a suspicious covariance, skipping", gene)
    #     status = CalculationStatus.SINGULAR_COVARIANCE
    #     return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status

    w = float(dot(dot(zscores, inv), zscores))
    chi2_p = stats.chi2.sf(w, n_indep)

    tmi = numpy.trace(numpy.dot(matrix,inv))

    # if we got to this point, we are  ok-ish. The chi distribution might have been unable to calculate the pvalue because it is too small...
    if chi2_p == 0:
        status = CalculationStatus.INSUFFICIENT_NUMERICAL_RESOLUTION
    else:
        status = CalculationStatus.OK

    pvalue = chi2_p

    if DEBUG:
        from IPython import embed; embed();

    return g, pvalue, n, n_indep, p_i_best, p_i_worst, eigen_max, eigen_min, z_min, z_max, z_mean, z_sd, tmi, status

def format_results(results):
    columns = ["gene", "pvalue", "n", "n_indep", "p_i_best", "p_i_worst", "eigen_max", "eigen_min", "z_min", "z_max", "z_mean", "z_sd", "tmi",  "status"]
    results = Utilities.to_dataframe(results, columns)
    results = results.sort_values(by=["pvalue", "status"])
    results = results.fillna("NA")
    if DEBUG:
        print "b"
        from IPython import embed; embed()
        exit()
    return results
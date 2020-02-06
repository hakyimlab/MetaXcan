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
    def get_gene_name(self, gene): raise  Exceptions.NotImplemented("Context: get_gene_name")
    def check(self): raise Exceptions.NotImplemented("Context: check")

class ContextMixin(object):
    def __init__(self):
        self.metaxcan_results_manager = None
        self.matrix_manager = None
        self.cutoff = None
        self.gene_names = None
        self.trimmed_ensemble_id = None

    def get_metaxcan_zscores(self, gene):
        if self.trimmed_ensemble_id and "." in gene:
            gene = gene.split(".")[0]
        results = self.metaxcan_results_manager.results_for_gene(gene)
        return results

    def get_model_matrix(self, gene, tissues):
        return self.matrix_manager.get(gene, tissues)

    def get_cutoff(self, matrix):
        return self.cutoff(matrix)

    def get_gene_name(self, gene):
        return self.gene_names[gene]

    def get_trimmed_ensemble_id(self):
        return self.trimmed_ensemble_id

    def _process_genes(self, genes):
        if self.trimmed_ensemble_id:
            g = {t.gene.split(".")[0]: t.gene_name for t in genes.itertuples()}
        else:
            g = {t.gene: t.gene_name for t in genes.itertuples()}
        return g

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

def joint_analysis(context, gene):
    g, g_n, pvalue, n, n_indep, p_i_best, t_i_best, p_i_worst, t_i_worst, eigen_max, eigen_min, eigen_min_kept, z_min, z_max, z_mean, z_sd, tmi, status \
        = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, CalculationStatus.NO_DATA
    g = gene.split(".")[0] if context.get_trimmed_ensemble_id() else gene

    g_n = context.get_gene_name(g)

    ####################################################################################################################
    zscores, tissue_labels = context.get_metaxcan_zscores(gene)
    if not zscores or len(zscores) == 0:
        status = CalculationStatus.NO_METAXCAN_RESULTS
        return g, g_n, pvalue, n, n_indep, p_i_best, t_i_best, p_i_worst, t_i_worst, eigen_max, eigen_min, eigen_min_kept, z_min, z_max, z_mean, z_sd, tmi, status
    n = len(zscores)
    z_min = numpy.min(zscores)
    z_max = numpy.max(zscores)
    z_mean = numpy.mean(zscores)
    if (len(zscores)>1):
        z_sd = numpy.std(zscores, ddof=1)

    ####################################################################################################################
    labels, matrix = context.get_model_matrix(gene, tissue_labels)

    if not labels or len(labels) == 0:
        status = CalculationStatus.NO_PRODUCT
        return g, g_n, pvalue, n, n_indep, p_i_best, t_i_best, p_i_worst, t_i_worst, eigen_max, eigen_min, eigen_min_kept, z_min, z_max, z_mean, z_sd, tmi, status

    # also, check that the matrix actually makes sense. We are currently returning it just in case but matrices with complex covariance are suspicious.
    e, v = numpy.linalg.eigh(matrix)
    if numpy.imag(e).any():
        status = CalculationStatus.COMPLEX_COVARIANCE
        e = numpy.real(e)
        eigen_max, eigen_min = numpy.max(e), numpy.min(e)
        return g, g_n, pvalue, n, n_indep, p_i_best, t_i_best, p_i_worst, t_i_worst, eigen_max, eigen_min, eigen_min_kept, z_min, z_max, z_mean, z_sd, tmi, status

    # If no eigenvalue satisfies our cutoff criteria, at least the first component will be used
    # Note there is a slight numerical mismatch between the resolution in eigh and the svd
    cutoff = context.get_cutoff(matrix)

    _d = {tissue_labels[i]:zscores[i] for i in range(0, len(tissue_labels))}
    zscores = array([_d[l] for l in labels])
    inv, n_indep, eigen = Math.capinv(matrix, cutoff, context.epsilon)

    eigen_max, eigen_min = numpy.max(eigen), numpy.min(eigen)
    eigen_min_kept = numpy.min([x for x in eigen[0:n_indep]])

    _absz = numpy.abs(zscores)
    _maxzi = numpy.argmax(_absz)
    max_z = _absz[_maxzi]
    p_i_best = 2*stats.norm.sf(max_z)
    t_i_best = labels[_maxzi]

    _minzi = numpy.argmin(_absz)
    min_z = _absz[_minzi]
    p_i_worst = 2*stats.norm.sf(min_z)
    t_i_worst = labels[_minzi]

    #TODO: implement a better heuristic
    try:
        eigen_w, eigen_v = numpy.linalg.eigh(inv)
    except:
        #WTCCC 'ENSG00000204560.5'
        logging.log(8, "Problems with inverse for %s, skipping", gene)
        status = CalculationStatus.INVERSE_ERROR
        return g, g_n, pvalue, n, n_indep, p_i_best, t_i_best, p_i_worst, t_i_worst, eigen_max, eigen_min, eigen_min_kept, z_min, z_max, z_mean, z_sd, tmi, status

    ####################################################################################################################
    w = float(dot(dot(zscores, inv), zscores))
    chi2_p = stats.chi2.sf(w, n_indep)

    tmi = numpy.trace(numpy.dot(matrix,inv))

    # if we got to this point, we are  ok-ish. The chi distribution might have been unable to calculate the pvalue because it is too small...
    if chi2_p == 0:
        status = CalculationStatus.INSUFFICIENT_NUMERICAL_RESOLUTION
    else:
        status = CalculationStatus.OK

    pvalue = chi2_p

    return g, g_n, pvalue, n, n_indep, p_i_best, t_i_best, p_i_worst, t_i_worst, eigen_max, eigen_min, eigen_min_kept, z_min, z_max, z_mean, z_sd, tmi, status

def format_results(results):
    columns = ["gene", "gene_name", "pvalue", "n", "n_indep", "p_i_best", "t_i_best", "p_i_worst", "t_i_worst", "eigen_max", "eigen_min", "eigen_min_kept", "z_min", "z_max", "z_mean", "z_sd", "tmi",  "status"]
    results = Utilities.to_dataframe(results, columns)
    results = results.sort_values(by=["pvalue", "status"])
    results = results.fillna("NA")
    return results
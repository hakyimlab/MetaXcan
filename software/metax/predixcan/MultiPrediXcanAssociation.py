import logging
import pandas

from patsy import dmatrices
import numpy
from numpy import dot as _dot, diag as _diag
import statsmodels.api as sm

from .. import Exceptions
from ..misc import Math

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract Multi Tissue PrediXcan context")
    def get_genes(self): raise  Exceptions.NotImplemented("MT Predixcan Context: get_genes")
    def expression_for_gene(self, gene): raise  Exceptions.NotImplemented("MT Predixcan Context: expression for genes")
    def get_pheno(self): raise  Exceptions.NotImplemented("MT Predixcan Context: get pheno")
    def get_mode(self): raise Exceptions.NotImplemented("MT Predixcan Context: get mode")
    def get_covariates(self): raise Exceptions.NotImplemented("MT Predixcan Context: get covariates")

class MTPF(object):
    """Multi-Tissue PrediXcan format"""
    GENE = 0
    PVALUE =1
    N_MODELS = 2
    N_SAMPLES = 3
    BEST_GWAS_P = 4
    BEST_GWAS_M = 5
    WORST_GWAS_P = 6
    WORST_GWAS_M = 7
    STATUS = 8
    N_USED=9
    MAX_EIGEN=10
    MIN_EIGEN=11
    MIN_EIGEN_KEPT=12

    K_GENE = "gene"
    K_PVALUE = "pvalue"
    K_N_MODELS = "n_models"
    K_N_SAMPLES = "n_samples"
    K_BEST_GWAS_P = "p_i_best"
    K_BEST_GWAS_M = "m_i_best"
    K_WORST_GWAS_P = "p_i_worst"
    K_WORST_GWAS_M = "m_i_worst"
    K_STATUS = "status"
    K_N_USED="n_used"
    K_MAX_EIGEN="max_eigen"
    K_MIN_EIGEN = "min_eigen"
    K_MIN_EIGEN_KEPT = "min_eigen_kept"

    order=[(GENE,K_GENE), (PVALUE, K_PVALUE), (N_MODELS,K_N_MODELS), (N_SAMPLES, K_N_SAMPLES), (BEST_GWAS_P, K_BEST_GWAS_P), (BEST_GWAS_M, K_BEST_GWAS_M), (WORST_GWAS_P, K_WORST_GWAS_P), (WORST_GWAS_M, K_WORST_GWAS_M),
           (STATUS, K_STATUS), (N_USED, K_N_USED), (MAX_EIGEN,K_MAX_EIGEN), (MIN_EIGEN, K_MIN_EIGEN), (MIN_EIGEN_KEPT,K_MIN_EIGEN_KEPT)]

########################################################################################################################
class MTPStatus(object):
    K_MLE_DID_NOT_CONVERGE="MLE_did_not_converge"

class MTPMode(object):
    K_LINEAR= "linear"
    K_LOGISTIC="logistic"

    K_MODES = [K_LINEAR, K_LOGISTIC]

def _ols_pvalue(result): return  result.f_pvalue
def _logit_pvalue(result): return result.llr_pvalue

def _ols_fit(model): return model.fit()
def _logit_fit(model): return model.fit(disp=False, maxiter=100)

def _ols_status(result): return None
def _logit_status(result):
    if not result.mle_retvals['converged']:
        return MTPStatus.K_MLE_DID_NOT_CONVERGE

    return None


K_METHOD="method"
K_FIT="fit"
K_PVALUE="pvalue"
K_STATUS="status"

_mode = {
    MTPMode.K_LINEAR:{
        K_METHOD:sm.OLS,
        K_FIT:_ols_fit,
        K_PVALUE:_ols_pvalue,
        K_STATUS:_ols_status
    },
    MTPMode.K_LOGISTIC:{
        K_METHOD: sm.Logit,
        K_FIT: _logit_fit,
        K_PVALUE: _logit_pvalue,
        K_STATUS: _logit_status
    }
}

########################################################################################################################

def _acquire(gene, context):
    e = context.expression_for_gene(gene)
    e = pandas.DataFrame(e)

    e["pheno"] = context.get_pheno()
    e_ = e.dropna()
    # discard columns where expression is zero for all individuals. Mostly happens in UKB diseases.
    e_ = e_[e_.columns[~(e_ == 0).all()]]

    model_keys = list(e_.columns.values)
    model_keys.remove("pheno")

    return model_keys, e_

def _design_matrices(e_, keys, context):
    formula = "pheno ~ {}".format(" + ".join(keys))
    y, X = dmatrices(formula, data=e_, return_type="dataframe")
    return y, X

def _pvalues(result, context):
    return  result.pvalues[result.pvalues.index[1:]]

def _get_pc_input(e_, model_keys):
    Xc = []
    _mk = []
    for key in model_keys:
        x = Math.standardize(e_[key])
        if x is not None:
            Xc.append(x)
            _mk.append(key)
    return Xc, _mk

def _pca_data(e_, model_keys, pc_filter):
    if e_.shape[1] == 2:
        return e_, model_keys, model_keys, 1, 1, 1, 1, 1
    #numpy.svd can't handle typical data size in UK Biobank. So we do PCA through the covariance matrix
    # That is: we compute ths SVD of a covariance matrix, and use those coefficients to get the SVD of input data
    # Shamelessly designed from https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
    # In numpy.cov, each row is a variable and each column an observation. Exactly opposite to standard PCA notation: it is transposed, then.
    Xc_t, original_keys = _get_pc_input(e_, model_keys)
    k = numpy.cov(Xc_t)
    u, s, vt = numpy.linalg.svd(k)
    # we want to keep only those components with significant variance, to reduce dimensionality
    selected = pc_filter(s)

    variance = s[selected]
    vt_projection = vt[selected]
    Xc_t_ = _dot(vt_projection, Xc_t)
    pca_keys = ["pc{}".format(i) for i in range(0, len(selected))]
    _data = {pca_keys[i]:x for i,x in enumerate(Xc_t_)}
    _data["pheno"] = e_.pheno
    pca_data = pandas.DataFrame(_data)
    return pca_data, pca_keys, original_keys, numpy.max(s), numpy.min(s), numpy.min(s[selected]), vt_projection, variance

def _coefs(result, vt_projection, model_keys):
    coefs = result.params[1:].to_frame().reset_index().rename(columns={"index": "variable", 0: "param"})
    if vt_projection is not None:
        v = numpy.dot(coefs.param.values,vt_projection)
        k = model_keys
        coefs = pandas.DataFrame({"variable": k, "param": v})
    return coefs

def multi_predixcan_association(gene_, context, callbacks=None):
    gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst, status, n_used, max_eigen, min_eigen, min_eigen_kept = None, None, None, None, None, None, None, None, None, None, None, None, None
    gene = gene_

    model_keys, e_ = _acquire(gene_, context)
    n_models = len(model_keys)
    pc_filter = context.get_pc_filter()
    try:
        if pc_filter is not None:
            e_, model_keys, original_models, max_eigen, min_eigen, min_eigen_kept, vt_projection, variance = _pca_data(e_, model_keys, pc_filter)
        else:
            original_models = model_keys
            vt_projection = None
            variance = None

        n_used = len(model_keys)
        y, X = _design_matrices(e_, model_keys, context)
        specifics =  _mode[context.get_mode()]
        model = specifics[K_METHOD](y, X)
        result = specifics[K_FIT](model)

        n_samples = e_.shape[0]

        p_i_ = _pvalues(result, context)
        p_i_best = p_i_.min()
        m_i_best = p_i_.idxmin()
        p_i_worst = p_i_.max()
        m_i_worst = p_i_.idxmax()

        pvalue = specifics[K_PVALUE](result)
        status = specifics[K_STATUS](result)

        if callbacks:
            coefs = _coefs(result, vt_projection, original_models)
            for callback in callbacks:
                callback(gene, model, result, vt_projection, variance, original_models, coefs)
    except Exception as ex:
        status = str(ex).replace(" ", "_").replace(",", "_")

    return gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst, status, n_used, max_eigen, min_eigen, min_eigen_kept

def dataframe_from_results(results, context):
    results = list(zip(*results))
    if len(results) == 0:
        return pandas.DataFrame({key:[] for order,key in MTPF.order})

    ORDER = MTPF.order if (context is not None and context.get_pc_filter() is not None) else MTPF.order[0:MTPF.N_USED]
    r = pandas.DataFrame({key: results[order] for order, key in ORDER})
    r = r[[key for order,key in ORDER]]
    return r

class SaveCoefs(object):
    def __init__(self):
        self.coefs = []

    def __call__(self, gene, model, result, vt_projection, variance, model_keys, coefs):
        self.coefs.append(coefs.assign(gene=gene))

    def get(self):
        return pandas.concat(self.coefs)

class SaveLoadings(object):
    def __init__(self):
        self.loadings = []

    def __call__(self, gene, model, result, vt_projection, variance,  model_keys, coefs):
        results=[]
        if type(vt_projection) == int or type(vt_projection) == float:
            results.append((gene, "pc0", model_keys[0], variance))
        else:
            for i in range(0, vt_projection.shape[0]):
                pc = "pc{}".format(i)
                l = numpy.sqrt(variance[i])*vt_projection[i]
                for j in range(0, l.shape[0]):
                    results.append((gene,pc,model_keys[j], l[j]))

        self.loadings.append(pandas.DataFrame(results, columns=["gene", "pc", "tissue", "weight"]))

    def get(self):
        return pandas.concat(self.loadings)

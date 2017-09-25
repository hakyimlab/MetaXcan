import logging
import pandas

from patsy import dmatrices
import statsmodels.api as sm

from .. import Exceptions

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

    K_GENE = "gene"
    K_PVALUE = "pvalue"
    K_N_MODELS = "n_models"
    K_N_SAMPLES = "n_samples"
    K_BEST_GWAS_P = "p_i_best"
    K_BEST_GWAS_M = "m_i_best"
    K_WORST_GWAS_P = "p_i_worst"
    K_WORST_GWAS_M = "m_i_worst"
    K_STATUS = "status"

    order=[(GENE,K_GENE), (PVALUE, K_PVALUE), (N_MODELS,K_N_MODELS), (N_SAMPLES, K_N_SAMPLES), (BEST_GWAS_P, K_BEST_GWAS_P), (BEST_GWAS_M, K_BEST_GWAS_M), (WORST_GWAS_P, K_WORST_GWAS_P), (WORST_GWAS_M, K_WORST_GWAS_M),
           (STATUS, K_STATUS),]

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
    if context.get_covariates() is not None:
        # If context used covariates, then the intercept was regressed out with them
        formula += " - 1"
    y, X = dmatrices(formula, data=e_, return_type="dataframe")
    return y, X

def _pvalues(result, context):
    if context.get_covariates() is None:
        return  result.pvalues[result.pvalues.index[1:]]

    return result.pvalues

def multi_predixcan_association(gene_, context):
    gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst, status = None, None, None, None, None, None, None, None, None
    gene = gene_

    model_keys, e_ = _acquire(gene_, context)
    n_models = len(model_keys)

    try:
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
    except Exception as ex:
        status = ex.message.replace(" ", "_").replace(",", "_")

    return gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst, status

def dataframe_from_results(results):
    results = zip(*results)
    if len(results) == 0:
        return pandas.DataFrame({key:[] for order,key in MTPF.order})

    r = pandas.DataFrame({key: results[order] for order, key in MTPF.order})
    r = r[[key for order,key in MTPF.order]]
    return r
import logging
import pandas

from patsy import dmatrices
import statsmodels.api as sm

from .. import Exceptions

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract Multi Tissue PrediXcan context")
    def get_genes(self): raise  Exceptions.NotImplemented("Predixcan Context: get_genes")
    def expression_for_gene(self, gene): raise  Exceptions.NotImplemented("Predixcan Context: expression for genes")
    def get_pheno(self): raise  Exceptions.NotImplemented("Predixcan Context: get pheno")
    def get_covariates(self): raise Exceptions.NotImplemented("Predixcan Context: get covariates")

class PF(object):
    """Multi-Tissue PrediXcan format"""
    GENE = 0
    EFFECT_SIZE =1
    SE=2
    ZSCORE=3
    PVALUE =4
    N_SAMPLES = 5
    STATUS = 6

    K_GENE = "gene"
    K_EFFECT_SIZE = "effect"
    K_SE = "se"
    K_ZSCORE = "zscore"
    K_PVALUE = "pvalue"
    K_N_SAMPLES = "n_samples"
    K_STATUS = "status"

    order=[(GENE,K_GENE), (EFFECT_SIZE,K_EFFECT_SIZE), (SE,K_SE), (ZSCORE, K_ZSCORE), (PVALUE, K_PVALUE),  (N_SAMPLES, K_N_SAMPLES),
           (STATUS, K_STATUS),]

########################################################################################################################
class PStatus(object):
    K_MLE_DID_NOT_CONVERGE="MLE_did_not_converge"

class PMode(object):
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
        return PStatus.K_MLE_DID_NOT_CONVERGE

    return None

K_METHOD="method"
K_FIT="fit"
K_PVALUE="pvalue"
K_STATUS="status"

_mode = {
    PMode.K_LINEAR:{
        K_METHOD:sm.OLS,
        K_FIT:_ols_fit,
        K_PVALUE:_ols_pvalue,
        K_STATUS:_ols_status
    },
    PMode.K_LOGISTIC:{
        K_METHOD: sm.Logit,
        K_FIT: _logit_fit,
        K_PVALUE: _logit_pvalue,
        K_STATUS: _logit_status
    }
}

########################################################################################################################

def _acquire(gene, context):
    e = context.expression_for_gene(gene)
    e = pandas.DataFrame({"expression":e})

    e["pheno"] = context.get_pheno()
    e_ = e.dropna()
    # discard columns where expression is zero for all individuals. Mostly happens in UKB diseases.
    e_ = e_[e_.columns[~(e_ == 0).all()]]

    model_keys = list(e_.columns.values)
    model_keys.remove("pheno")

    return model_keys, e_

def _design_matrices(e_,  context):
    formula = "pheno ~ expression"
    y, X = dmatrices(formula, data=e_, return_type="dataframe")
    return y, X

def _results(result, context):
    idx = 1

    effect = result.params[idx]
    se = result.bse[idx]
    zscore = result.tvalues[idx]
    pvalue = result.pvalues[idx]

    return effect, se, zscore, pvalue

def predixcan_association(gene_, context):
    gene, effect_size, se, zscore, pvalue, n_samples, status = None, None, None, None, None, None, None
    gene = gene_
    model_keys, e_ = _acquire(gene_, context)

    try:
        y, X = _design_matrices(e_,  context)
        specifics =  _mode[context.get_mode()]
        model = specifics[K_METHOD](y, X)
        result = specifics[K_FIT](model)
        n_samples = e_.shape[0]
        effect_size, se, zscore, pvalue = _results(result, context)
        status = specifics[K_STATUS](result)
    except Exception as ex:
        status = str(ex).replace(" ", "_").replace(",", "_")

    return gene, effect_size, se, zscore, pvalue, n_samples, status

def dataframe_from_results(results):
    results = list(zip(*results))
    if len(results) == 0:
        return pandas.DataFrame({key:[] for order,key in PF.order})

    r = pandas.DataFrame({key: results[order] for order, key in PF.order})
    r = r[[key for order,key in PF.order]]
    return r
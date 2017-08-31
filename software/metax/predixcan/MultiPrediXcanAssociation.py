from .. import Exceptions

from patsy import dmatrices
import statsmodels.api as sm

import pandas

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract Multi Tissue PrediXcan context")
    def get_genes(self): raise  Exceptions.NotImplemented("MT Predixcan Context: get_genes")
    def expression_for_gene(self, gene): raise  Exceptions.NotImplemented("MT Predixcan Context: expression for genes")
    def get_pheno(self): raise  Exceptions.NotImplemented("MT Predixcan Context: get pheno")
    def get_mode(self): raise Exceptions.NotImplemented("MT Predixcan Context: get mode")

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

    K_GENE = "gene"
    K_PVALUE = "pvalue"
    K_N_MODELS = "n_models"
    K_N_SAMPLES = "n_samples"
    K_BEST_GWAS_P = "p_i_best"
    K_BEST_GWAS_M = "m_i_best"
    K_WORST_GWAS_P = "p_i_worst"
    K_WORST_GWAS_M = "m_i_worst"

    order=[(GENE,K_GENE), (PVALUE, K_PVALUE), (N_MODELS,K_N_MODELS), (N_SAMPLES, K_N_SAMPLES), (BEST_GWAS_P, K_BEST_GWAS_P), (BEST_GWAS_M, K_BEST_GWAS_M), (WORST_GWAS_P, K_WORST_GWAS_P), (WORST_GWAS_M, K_WORST_GWAS_M),]

########################################################################################################################
class MTPMode(object):
    K_LINEAL="lineal"
    K_LOGISTIC="logistic"

    K_MODES = [K_LINEAL, K_LOGISTIC]

def _ols_pvalue(result): return  result.f_pvalue
def _logit_pvalue(result): return result.llr_pvalue

def _ols_fit(model): return model.fit()
def _logit_fit(model): return model.fit(disp=False)

K_METHOD="method"
K_FIT="fit"
K_PVALUE="pvalue"

_mode = {
    MTPMode.K_LINEAL: {
        K_METHOD:sm.OLS,
        K_FIT:_ols_fit,
        K_PVALUE:_ols_pvalue
    },
    MTPMode.K_LOGISTIC:{
        K_METHOD:sm.Logit,
        K_FIT:_logit_fit,
        K_PVALUE:_logit_pvalue
    }
}

########################################################################################################################

def multi_predixcan_association(gene_, context):
    gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst = None, None, None, None, None, None, None, None
    gene = gene_

    e = context.expression_for_gene(gene)
    n_models = len(e.keys())
    model_keys = sorted(e.keys())

    e = pandas.DataFrame(e)
    e["y"] = context.get_pheno()
    e_ = e.dropna()

    y, X = dmatrices("y ~ {}".format(" + ".join(model_keys)), data=e_, return_type="dataframe")
    specifics =  _mode[context.get_mode()]
    model = specifics[K_METHOD](y, X)
    result = specifics[K_FIT](model)

    n_samples = e_.shape[0]

    p_i_ = result.pvalues[result.pvalues.index[1:]]
    p_i_best = p_i_.min()
    m_i_best = p_i_.idxmin()
    p_i_worst = p_i_.max()
    m_i_worst = p_i_.idxmax()

    pvalue = specifics[K_PVALUE](result)

    return gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst

def dataframe_from_results(results):
    results = zip(*results)
    if len(results) == 0:
        return pandas.DataFrame({key:[] for order,key in MTPF.order})

    r = pandas.DataFrame({key: results[order] for order, key in MTPF.order})
    r = r[[key for order,key in MTPF.order]]
    return r
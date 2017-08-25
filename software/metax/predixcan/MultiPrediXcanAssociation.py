from .. import Exceptions

from patsy import dmatrices
import statsmodels.api as sm

import pandas

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract Multi Tissue PrediXcan context")
    def get_genes(self): raise  Exceptions.NotImplemented("MT Predixcan Context: get_genes")
    def expression_for_gene(self, gene): raise  Exceptions.NotImplemented("MT Predixcan Context: expression for genes")
    def get_pheno(self): raise  Exceptions.NotImplemented("MT Predixcan Context: get pheno")

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


def multi_predixcan_association(gene_, context):
    gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst = None, None, None, None, None, None, None, None
    gene = gene_

    e = context.expression_for_gene(gene)
    n_models = len(e.keys())
    model_keys = sorted(e.keys())

    #{specific
    e = pandas.DataFrame(e)
    e["y"] = context.get_pheno()
    e_ = e.dropna()
    y, X = dmatrices("y ~ {}".format(" + ".join(model_keys)), data=e_, return_type="dataframe")
    model = sm.OLS(y, X)
    result = model.fit()

    n_samples = e_.shape[0]

    p_i_ = result.pvalues[result.pvalues.index[1:]]
    p_i_best = p_i_.min()
    m_i_best = p_i_.idxmin()
    p_i_worst = p_i_.max()
    m_i_worst = p_i_.idxmax()

    pvalue = result.f_pvalue
    #specific}

    return gene, pvalue, n_models, n_samples, p_i_best, m_i_best, p_i_worst,  m_i_worst

def dataframe_from_results(results):
    results = zip(*results)
    if len(results) == 0:
        return pandas.DataFrame({key:[] for order,key in MTPF.order})

    r = pandas.DataFrame({key: results[order] for order, key in MTPF.order})
    r = r[[key for order,key in MTPF.order]]
    return r
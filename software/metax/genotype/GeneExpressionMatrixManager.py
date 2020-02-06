import numpy
from numpy import dot as _d

from .. import MatrixManager2
from .. import MatrixManager

###############################################################################
class GeneExpressionMatrixManager(object):
    def __init__(self, snp_covariance, model_manager):
        self.snp_covariance_manager = MatrixManager2.MatrixManager2(snp_covariance, MatrixManager.GENE_SNP_COVARIANCE_DEFINITION)
        self.model_manager = model_manager

    def get(self, gene, tissues):
        models = self.model_manager.get_models(gene)
        models = models.loc[tissues]
        tissues = sorted(set(models.index.get_level_values(0).values))
        variances = _get_variances(models, self.snp_covariance_manager, gene)
        tissues, matrix = _build_matrix(gene, models, self.snp_covariance_manager, tissues, variances, _get_coef)
        return tissues, matrix

def _build_matrix(gene, models, matrix_manager, tissues, variances, coef_method):
    coefs, _tissues = _build_matrix_entries(gene, models, matrix_manager, tissues, variances, coef_method)
    matrix = MatrixManager._to_matrix(coefs, _tissues)
    return _tissues, matrix

def _build_matrix_entries(gene, models, matrix_manager, tissues, variances, coef_method):
    coefs = {}
    _t = set()
    _tissues = []
    for i in range(0, len(tissues)):
        for j in range(i, len(tissues)):
            t1 = tissues[i]
            t2 = tissues[j]
            if not t1 in coefs: coefs[t1] = {}
            if not t2 in coefs: coefs[t2] = {}

            value = coef_method(gene, models, matrix_manager, variances, t1, t2)

            if value is None or value is numpy.nan:
                continue

            coefs[t1][t2] = value
            coefs[t2][t1] = value
            if not t1 in _t:
                _t.add(t1)
                _tissues.append(t1)
    return coefs, _tissues

def _get_variances(models, matrix_manager, gene):
    tissues = models.index.get_level_values(0).values
    variances = {t:_get_variance(models, matrix_manager, gene, t) for t in tissues}
    variances = {k:v for k,v in variances.items() if v is not None}
    return variances

def _get_variance(models, matrix_manager, gene, tissue):
    model = models.loc[tissue]
    snps = set(model.index.get_level_values(0).values)

    #Remember that only those snps with data in the GWAS will get loaded.
    snps, matrix = matrix_manager.get(gene, snps, strict_whitelist=False)
    if len(snps) == 0:
        return None

    weights = model.loc[snps].weight.values
    variance = _d(_d(weights, matrix), weights)
    variance = numpy.float64(variance)
    return variance

def _get_coef(gene, models,  matrix_manager, variances, t1, t2):
    model_1 = models.loc[t1]
    snps_1 = set(model_1.index.get_level_values(0))
    if len(snps_1) == 0:
        return None

    model_2 = models.loc[t2]
    snps_2 = set(model_2.index.get_level_values(0))
    if len(snps_2) == 0:
        return None

    s1, s2, matrix = matrix_manager.get_2(gene, snps_1, snps_2)

    if len(s1) == 0 or len(s2) == 0:
        return None

    w1 = model_1.loc[s1].weight.values
    w2 = model_2.loc[s2].weight.values

    denom = numpy.float64(numpy.sqrt(variances[t1] * variances[t2]))
    if denom == 0:
        return numpy.nan

    num = numpy.float64(_d(_d(w1, matrix), w2))
    coef = num/denom
    return coef

###############################################################################
# Like the above but faster and dirtier
class _GeneExpressionMatrixManager(object):
    """Like 'GeneExpressionMatrixManager', but relying on faster implementation."""
    def __init__(self, snp_covariance, model_manager):
        self.snp_covariance_manager = MatrixManager2.MatrixManager2(snp_covariance, MatrixManager.GENE_SNP_COVARIANCE_DEFINITION)
        self.model_manager = model_manager

    def get(self, gene, tissues):
        t = set(tissues)
        models = self.model_manager.get_models(gene)
        models = {k:v for k,v in models.items() if k in t}
        tissues = sorted(set(models.keys()))
        variances = _get_variances_2(models, self.snp_covariance_manager, gene)
        tissues, matrix = _build_matrix(gene, models, self.snp_covariance_manager, tissues, variances, _get_coef_2)
        return tissues, matrix

def _get_variances_2(models, matrix_manager, gene):
    tissues = list(models.keys())
    variances = {t:_get_variance_2(models, matrix_manager, gene, t) for t in tissues}
    variances = {k:v for k,v in variances.items() if v is not None}
    return variances

def _get_variance_2(models, matrix_manager, gene, tissue):
    model = models[tissue]
    snps = set(model.keys())

    #Remember that only those snps with data in the GWAS will get loaded.
    snps, matrix = matrix_manager.get(gene, snps, strict_whitelist=False)
    if len(snps) == 0:
        return None

    weights = numpy.array([model[x] for x in snps], dtype=numpy.float64)
    variance = _d(_d(weights, matrix), weights)
    variance = numpy.float64(variance)
    return variance

def _get_coef_2(gene, models,  matrix_manager, variances, t1, t2):
    model_1 = models[t1]
    snps_1 = set(model_1.keys())
    if len(snps_1) == 0:
        return None

    model_2 = models[t2]
    snps_2 = set(model_2.keys())
    if len(snps_2) == 0:
        return None

    s1, s2, matrix = matrix_manager.get_2(gene, snps_1, snps_2)

    if len(s1) == 0 or len(s2) == 0:
        return None

    w1 = numpy.array([model_1[x] for x in s1], dtype=numpy.float64)
    w2 = numpy.array([model_2[x] for x in s2], dtype=numpy.float64)

    denom = numpy.float64(numpy.sqrt(variances[t1] * variances[t2]))
    if denom == 0:
        return numpy.nan

    num = numpy.float64(_d(_d(w1, matrix), w2))
    coef = num/denom
    return coef

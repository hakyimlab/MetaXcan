import numpy
import logging

from numpy import dot as _d

from .. import MatrixManager
from .. import MatrixManager2
from .. import Utilities
from ..misc import Math


class GenotypeAnalysisContext(object):
    def __init__(self, metadata, dosage, model_manager, standardize=False):
        self.metadata = metadata
        self.dosage = numpify(dosage)
        self.model_manager = model_manager
        if standardize: self.standardise_data()

    def standardise_data(self):
        d = {}
        rejected = set()
        for k,v in self.dosage.iteritems():
            v = Math.standardize(v)
            if v is None:
                rejected.add(k)
                continue
            d[k] = v
        self.dosage = d
        self.metadata = self.metadata[~self.metadata.rsid.isin(rejected)]

    def get_genes(self):
        logging.log(7,"getting snps in dosage")
        dosage_snps = set(self.dosage.keys())
        logging.log(7, "getting genes for snps")
        genes = self.model_manager.get_implicated_genes(dosage_snps)
        logging.log(7, "returning genes")
        return genes

    def get_rsids(self, gene=None):
        dosage_snps = set(self.dosage.keys())
        model_snps = self.model_manager.get_rsids(gene)
        snps = dosage_snps.intersection(model_snps)
        return set(snps)

    def get_model_labels(self, gene=None):
        return self.model_manager.get_model_labels(gene)

    def get_models(self, gene=None):
        return self.model_manager.get_models(gene)

    def get_dosage(self, rsid):
        return self.dosage[rsid]

def numpify(d):
    return {k:numpy.array(v, dtype=numpy.float64) for k, v in d.iteritems()}


def get_prediction_variance(context, gene):
    labels = context.get_model_labels(gene)
    models = context.get_models(gene)
    snps = context.get_rsids(gene)

    results = []
    for label in labels:
        model = models.loc[label]
        model_snps = [x for x in model.index.values if x in snps]
        n_snps = len(model_snps)
        T = [context.get_dosage(l) * model.loc[l].weight for l in model_snps] #(X_0 * w_0, ... ,X_p * w_p)
        T = numpy.sum(T, axis=0) #Sum_l (X_l * w_l)
        # Mind the degrees of freedom! This way we will match  w_ * GAMMA * w, because covariance was calculated with ddof=1
        v = numpy.var(T, ddof=1)
        results.append((gene, label, v, n_snps))

    return results

VARIANCE_COLUMNS=["gene", "model", "variance", "n_snps"]
def format_prediction_variance_results(results):
    results = Utilities.to_dataframe(results, VARIANCE_COLUMNS)
    results = results.sort_values(by=["gene", "model"])
    results = results.fillna("NA")
    return results

def get_prediction_covariance(context, gene):
    snps = sorted(context.get_rsids(gene))
    X = [context.get_dosage(x) for x in snps]
    cov = numpy.cov(X)
    return gene, snps, cov

COVARIANCE_COLUMNS=["GENE", "RSID1", "RSID2", "VALUE"]
def format_prediction_covariance_results(results):
    flat = []
    for result in results:
        data = MatrixManager._flatten_matrix_data([result])
        flat.extend(data)

    data = Utilities.to_dataframe(flat, COVARIANCE_COLUMNS)
    data = data.fillna("NA")
    data = data.sort_values(by=COVARIANCE_COLUMNS[0:2])
    return data

class GeneExpressionMatrixManager(object):
    def __init__(self, snp_covariance, model_manager):
        self.snp_covariance_manager = MatrixManager2.MatrixManager2(snp_covariance, MatrixManager.GENE_SNP_COVARIANCE_DEFINITION)
        self.model_manager = model_manager

    def get(self, gene, tissues):
        models = self.model_manager.get_models(gene)
        models = models.loc[tissues]
        tissues, matrix = _build_matrix(gene, models, self.snp_covariance_manager)
        return tissues, matrix

def _build_matrix(gene, models, matrix_manager):
    tissues = sorted(set(models.index.get_level_values(0).values))
    variances = _get_variances(models, matrix_manager, gene)
    coefs = {}
    _t = set()
    _tissues = []
    for i in xrange(0, len(tissues)):
        for j in xrange(i, len(tissues)):
            t1 = tissues[i]
            t2 = tissues[j]
            if not t1 in coefs: coefs[t1] = {}
            if not t2 in coefs: coefs[t2] = {}

            value = _get_coef(gene, models, matrix_manager, variances, t1, t2)

            if not value:
                continue

            coefs[t1][t2] = value
            coefs[t2][t1] = value
            if not t1 in _t:
                _t.add(t1)
                _tissues.append(t1)

    matrix = MatrixManager._to_matrix(coefs, _tissues)
    return _tissues, matrix

def _get_variances(models, matrix_manager, gene):
    tissues = models.index.get_level_values(0).values
    variances = {t:_get_variance(models, matrix_manager, gene, t) for t in tissues}
    variances = {k:v for k,v in variances.iteritems() if v is not None}
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

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
        for k,v in self.dosage.items():
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
    return {k:numpy.array(v, dtype=numpy.float64) for k, v in d.items()}


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

import os
import pandas
import logging
from numpy import dot, transpose, cov

from .. import MatrixManager
from .. import Utilities
from . import Math

K_GENE_MODEL_HEADER=["gene", "model1", "model2", "value"]

class FeatureMatrixManager(object):
    def __init__(self, data, standardize=True):
        self.data = _build_data(data, standardize)
        self.columns = _get_columns(data)

    def get_feature_product(self, feature_key, center=False):
        X, labels = self._x(feature_key)
        product = dot(X, transpose(X)) if not center else cov(X) #mind, works the other way around to our notation
        return product, labels

    def get_feature_cov(self, feature_key):
        X, labels = self._x(feature_key)
        product = cov(X)
        return product, labels

    def _x(self, feature_key):
        vectors = self.data[feature_key]
        labels = sorted(vectors.keys())
        X = [vectors[m] for m in labels]
        return X, labels

    def save_covariances(self, path, column_names=K_GENE_MODEL_HEADER):
        data = []
        for key in self.columns:
            matrix, labels = self.get_feature_cov(key)
            flat = MatrixManager._flatten_matrix_data([(key, labels, matrix)])
            data.extend(flat)

        data = Utilities.to_dataframe(data, column_names)
        data = data.sort_values(by=["gene", "model1", "model2"])
        compression = "gzip" if "gz" in path else None
        data.to_csv(path, index=False, sep="\t", compression=compression)

#TODO: not necessarily "genes", but genetic features
def build_manager(folder, filters=["TW_*"], standardize=True):
    files = Utilities.target_files(folder, file_filters=filters)
    expressions = _load_features(files, _parse_expression_name)
    manager = FeatureMatrixManager(expressions, standardize)
    return manager

def _parse_expression_name(file):
    n = os.path.basename(os.path.normpath(file))
    n = n.replace("TW_","").replace(".expr.txt", "").replace("_0.5", "")
    return n

def _load_features(files, parse_func):
    results = {}
    for file in sorted(files):
        logging.log(9, "Loading feature file %s", file)
        e = pandas.read_csv(file, sep="\t")
        tag = parse_func(file)
        results[tag] = e
    return results

def _build_data(data, standardize):
    logging.log(9,"Building data")
    result = {}
    for k, df in data.iteritems():
        logging.log(9, "Processing %s", k)
        columns = df.columns.values
        for column in columns:
            if not column in result:
                result[column] = {}
            c = result[column]
            values =  df[column].values
            if standardize:
                values = Math.standardize(values)
            c[k] = values
    return result

def _get_columns(data):
    logging.info("Getting columns")
    columns = {x for k in data.keys() for x in data[k].columns}
    return columns
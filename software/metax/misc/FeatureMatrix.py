import os
import pandas
import logging
from numpy import dot, transpose

from .. import MatrixManager
from .. import Utilities

K_GENE_MODEL_HEADER=["gene", "model1", "model2", "value"]

class FeatureMatrixManager(object):
    def __init__(self, data):
        self.data = _build_data(data)
        self.columns = _get_columns(data)

    def get_feature_matrix(self, feature_key):
        vectors = self.data[feature_key]
        labels = sorted(vectors.keys())
        X = [vectors[m] for m in labels]
        product = dot(X, transpose(X)) # this is the other way around to our notation
        return product, labels

    def save(self, path, column_names=K_GENE_MODEL_HEADER):
        data = []
        for key in self.columns:
            matrix, labels = self.get_feature_matrix(key)
            flat = MatrixManager._flatten_matrix_data([(key, labels, matrix)])
            data.extend(flat)

        data = Utilities.to_dataframe(data, column_names)
        data = data.sort_values(by=["gene", "model1", "model2"])
        data.to_csv(path, index=False, sep="\t")

#TODO: not necessarily "genes", but genetic features
def build_manager(folder, filters=["TW_*"]):
    files = Utilities.target_files(folder, file_filters=filters)
    expressions = _load_features(files, _parse_expression_name)
    manager = FeatureMatrixManager(expressions)
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

def _build_data(data):
    logging.log(9,"Building data")
    result = {}
    for k, df in data.iteritems():
        logging.log(9, "Processing %s", k)
        columns = df.columns.values
        for column in columns:
            if not column in result:
                result[column] = {}
            c = result[column]
            c[k] = df[column].values
    return result

def _get_columns(data):
    logging.info("Getting columns")
    columns = {x for k in data.keys() for x in data[k].columns}
    return columns
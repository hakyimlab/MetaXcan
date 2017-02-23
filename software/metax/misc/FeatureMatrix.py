import os
import pandas
from numpy import dot, transpose
from numpy.linalg import inv

from .. import Utilities

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

    def save(self, path):
        pass

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
    files = files[0:2]
    results = {}
    for file in files:
        e = pandas.read_csv(file, sep="\t")
        tag = parse_func(file)
        results[tag] = e
    return results

def _build_data(data):
    result = {}
    for k, df in data.iteritems():
        columns = df.columns.values
        for column in columns:
            if not column in result:
                result[column] = {}
            c = result[column]
            c[k] = df[column]
    return result

def _get_columns(data):
    columns = {x for k in data.keys() for x in data[k].columns}
    return columns
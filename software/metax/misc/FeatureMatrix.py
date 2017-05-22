import pandas
import logging
import re
import gzip
from numpy import dot, transpose, cov

from .. import MatrixManager
from .. import Utilities
from . import Math

K_GENE_MODEL_HEADER=["gene", "model1", "model2", "value"]

class FeatureMatrixManager(object):
    def __init__(self, data, standardize=True):
        self.data = _build_data(data, standardize)
        self.columns = _get_columns(data)
        self.standardize = True

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

    def save_covariances(self, path, column_names=K_GENE_MODEL_HEADER, append=False):
        data = []
        for key in self.columns:
            matrix, labels = self.get_feature_cov(key)
            flat = MatrixManager._flatten_matrix_data([(key, labels, matrix)])
            data.extend(flat)

        data = Utilities.to_dataframe(data, column_names)
        data = data.sort_values(by=column_names[0:2])
        compression = "gzip" if "gz" in path else None
        if append:
            _open = gzip.open if compression == "gzip" else open
            with _open(path, 'a') as f:
                data.to_csv(f, header=False, index=False, sep="\t", compression=compression)
        else:
            data.to_csv(path, index=False, sep="\t", compression=compression)

#TODO: not necessarily "genes", but genetic features
def build_manager(folder, filters=["TW_*"], standardize=True, subset=None):
    files = Utilities.target_files(folder, file_filters=filters)
    expressions = _load_features(files, _parse_expression_name, subset)
    manager = FeatureMatrixManager(expressions, standardize)
    return manager

def _features(file):
    r = None
    with open(file) as f:
        r = f.readline().strip().split()
    return r

def features_in_folder(folder, filters=["TW_*"]):
    features = set()
    files = Utilities.target_files(folder, file_filters=filters)
    for file in files:
        l = _features(file)
        features.update(l)
    return features

_regexp = re.compile(".*TW_(.*)_0.5.expr.txt$")
def _parse_expression_name(file):
    return _regexp.match(file).group(1)

def _load_features(files, parse_func, subset):
    subset = set(subset) if subset else None
    results = {}
    for file in sorted(files):
        logging.log(9, "Loading feature file %s", file)
        features = None
        if subset:
            features = _features(file)
            features = [x for x in features if x in subset]

            if len(features) == 0:
                logging.info("No selected features for %s, skipping", file)
                continue

        e = pandas.read_csv(file, sep="\t", usecols=features)
        tag = parse_func(file)
        results[tag] = e
    return results

# dictionary[gene][tissue]
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
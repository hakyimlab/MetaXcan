"""
This is an alternative version of MatrixManager. This is way faster but consumes more memory.
The ordinary MatrixManager performs well enough when a gene's data is queried only once as in MetaXcan,
and is memory efficient. This one, on the contrary, organizes data in a way that is fast to query but takes up more memory
"""
import numpy
from . import MatrixManager

class MatrixManager2(MatrixManager.MatrixManagerBase):
    """
    Needs a dictionary mapping the header names to the keys ["model", "id1", "id2", "value"],
    """
    def __init__(self, d, definition):
        self.definition = definition
        MatrixManager._validate(d, definition)
        self.data = _build_data(d, definition)


    def get(self, key, whitelist=None, strict_whitelist=True):
        return _get(self.data, key, whitelist, strict_whitelist)

    def get_2(self, key, snps_1, snps_2):
        return _get_2(self.data, key, snps_1, snps_2)

    def model_labels(self):
        return set(self.data.keys())

    def n_ids(self, gene):
        if not gene in self.data:
            return numpy.nan
        d = self.data[gene]
        return len(list(d.keys()))

def _build_data(d, definition):
    #First element in a pandas tuple is the index.
    MODEL = d.columns.get_loc(definition[MatrixManager.K_MODEL]) + 1
    ID1 = d.columns.get_loc(definition[MatrixManager.K_ID1]) + 1
    ID2 = d.columns.get_loc(definition[MatrixManager.K_ID2]) + 1
    VALUE = d.columns.get_loc(definition[MatrixManager.K_VALUE]) + 1

    d = d.fillna("NA")
    r = {}
    for t in d.itertuples():
        model = t[MODEL]
        if not model in r: r[model] = {}

        id1 = t[ID1]
        id2 = t[ID2]
        value = t[VALUE]
        if value == "NA":
            continue

        m = r[model]
        if not id1 in m: m[id1] = {}
        if not id2 in m: m[id2] = {}
        m[id1][id2] = value
        m[id2][id1] = value
    return r

def _get(data, key, whitelist=None, strict_whitelist=True):
    if not key in data:
        return None,None

    d = data[key]
    ids = sorted(d.keys())

    if strict_whitelist and whitelist:
        MatrixManager._check_strict(whitelist, set(ids), key)

    if whitelist:
        ids = [x for x in ids if x in whitelist]

    covariance_matrix = MatrixManager._to_matrix(d, ids)
    return ids, covariance_matrix

def _get_2(data, key, id_1, id_2):
    if not key in data:
        return None,None

    d = data[key]
    snps = sorted(d.keys())
    _snps = set(snps)

    is1 = sorted([x for x in id_1 if x in _snps])
    is2 = sorted([x for x in id_2 if x in _snps])
    matrix = MatrixManager._to_matrix(d, is1, is2)
    return is1, is2, matrix
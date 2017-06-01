import pandas
import numpy
import Exceptions

K_MODEL="model"
K_ID1="id1"
K_ID2="id2"
K_VALUE="value"

GENE_SNP_COVARIANCE_DEFINITION = {
    K_MODEL:"GENE",
    K_ID1:"RSID1",
    K_ID2:"RSID2",
    K_VALUE:"VALUE"
}

class CDTF(object):
    """How we organize in memory data at the matrix"""
    MODEL=0
    ID1=1
    ID2=2
    VALUE=3

def load_matrix_manager(path, definition=GENE_SNP_COVARIANCE_DEFINITION):
    d = pandas.read_table(path, sep="\s+")
    m = MatrixManager(d, definition)
    return m

class MatrixManager(object):
    """
    Needs a dictionary mapping the header names to the keys ["model", "id1", "id2", "value"],
    """
    def __init__(self, d, definition):
        self.definition = definition
        _validate(d, definition)
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
        snps = self.data[gene]
        snps = _non_na(snps)
        snps = {x[CDTF.ID1] for x in snps}
        return len(snps)

def _validate(d, definition):
    MODEL_KEY = definition[K_MODEL]
    processed_genes = set()
    last_gene = None
    genes = d[MODEL_KEY]
    for g in genes:
        if g != last_gene:
            if g in processed_genes:
                msg = "Snp Matrix Entries for genes must be contiguous but %s was found in two different places" % (g)
                raise Exceptions.InvalidInputFormat(msg)
            processed_genes.add(g)
            last_gene = g

    if numpy.any(d.duplicated()):
        msg = "Duplicated SNP entries found"
        raise Exceptions.InvalidInputFormat(msg)

def _build_data(d, definition):
    MODEL_KEY=definition[K_MODEL]
    ID1_KEY=definition[K_ID1]
    ID2_KEY = definition[K_ID2]
    VALUE_KEY = definition[K_VALUE]

    d = d.fillna("NA")
    d[MODEL_KEY] = pandas.Categorical(d[MODEL_KEY], d[MODEL_KEY].drop_duplicates())  # speed things up!
    d = zip(d[MODEL_KEY].values, d[ID1_KEY].values, d[ID2_KEY].values, d[VALUE_KEY].values)
    r = {}
    for t in d:
        model = t[0]
        if not model in r:
            r[model] = []
        r[model].append(t)
    return r

def _get(data, key, whitelist=None, strict_whitelist=True):
    if not key in data:
        return None,None

    d = data[key]

    if strict_whitelist and whitelist:
        g, r1, r2, v = zip(*d)
        _check_strict(whitelist, set(r1), key)
        _check_strict(whitelist, set(r2), key)

    snps, entries = _rows_to_entries(d, key, whitelist)
    covariance_matrix = _to_matrix(entries, snps)
    return snps, covariance_matrix

def _rows_to_entries(d, key, whitelist):
    entries = {}
    _i = set()
    ids = []
    for row in d:
        id1 = row[CDTF.ID1]
        id2 = row[CDTF.ID2]
        if whitelist:
            if not id1 in whitelist: continue
            if not id2 in whitelist: continue

        value = row[CDTF.VALUE]
        if value == "NA": continue
        _check_value(value, key, id1, id2)

        if not id1 in entries: entries[id1] = {}
        if not id2 in entries: entries[id2] = {}

        entries[id1][id2] = value
        entries[id2][id1] = value

        if not id1 in _i:
            _i.add(id1)
            ids.append(id1)
    return ids, entries

def _check_strict(whitelist, snps, key):
    extra = {x for x in whitelist if not x in snps}
    if len(extra):
        msg = "SNPs in whitelist not in matrix for %s:%s" % (key, extra)
        raise Exceptions.InvalidArguments(msg)

def _check_value(value, key, id1, id2):
    try:
        float(value)
    except:
        msg = "Invalid value:{} for ({},{},{})".format(value, key, id1, id2)
        raise Exceptions.InvalidInputFormat(msg)

def _get_2(data, key, id_1, id_2):
    if not key in data:
        return None,None

    d = data[key]
    whitelist = set(id_1)|set(id_2)
    snps, entries = _rows_to_entries(d, key, whitelist)

    is1 = sorted([x for x in id_1 if x in snps])
    is2 = sorted([x for x in id_2 if x in snps])
    matrix = _to_matrix(entries, is1, is2)
    return is1, is2, matrix

def _to_matrix(entries, keys_i, keys_j=None):
    if not keys_j: keys_j = keys_i
    rows = []
    for key_i in keys_i:
        row = []
        rows.append(row)
        for key_j in keys_j:
            row.append(entries[key_i][key_j])
    matrix = numpy.matrix(rows, dtype=numpy.float64)
    return matrix

def _non_na(snps_data):
    return [x for x in snps_data if  x[CDTF.VALUE] != "NA"]

def _flatten_matrix_data(data):
    """data is expected to be a list of (name, id_labels, matrix) tuples"""
    results = []
    for name, id_labels, matrix in data:
        if len(id_labels) == 1:
            id = id_labels[0]
            results.append((name, id, id, float(matrix)))
            continue

        for i in xrange(0, len(id_labels)):
            for j in xrange(i, len(id_labels)):
                value = matrix[i][j]
                id1 = id_labels[i]
                id2 = id_labels[j]
                results.append((name, id1, id2, value))
    return results
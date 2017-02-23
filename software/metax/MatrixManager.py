import pandas
import numpy
import Exceptions

def load_matrix_manager(path):
    d = pandas.read_table(path, sep="\s+")
    m = MatrixManager(d)
    return m


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

class MatrixManager(object):
    """
    Needs a dictionary mapping the header names to the keys ["model", "id1", "id2", "value"],
    """
    def __init__(self, d, definition=GENE_SNP_COVARIANCE_DEFINITION):
        self.definition = definition
        _validate(d, definition)
        self.data = _build_data(d, definition)

    def get(self, key, whitelist=None, strict=True):
        return _get(self.data, key, whitelist, strict)

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
    d.GENE = pandas.Categorical(d.GENE, d.GENE.drop_duplicates())  # speed things up!
    d = zip(d[MODEL_KEY].values, d[ID1_KEY].values, d[ID2_KEY].values, d[VALUE_KEY].values)
    r = {}
    for t in d:
        model = t[0]
        if not model in r:
            r[model] = []
        r[model].append(t)
    return r

def _get(d, key, whitelist=None, strict=True):
    if not key in d:
        return None,None

    d = d[key]

    if whitelist is not None:
        g, r1, r2, v = zip(*d)
        snps = set(r1)
        whitelist = set(whitelist)
        if strict:
            extra = {x for x in whitelist if not x in snps}
            if len(extra):
                msg = "SNPs in whitelist not in matrix for %s:%s"%(key, extra)
                raise Exceptions.InvalidArguments(msg)
        d = [x for x in d if x[CDTF.ID1] in whitelist]

    _s = set()
    snps = []
    entries = {}
    for row in d:
        id1 = row[CDTF.ID1]
        id2 = row[CDTF.ID2]
        if not id1 in entries: entries[id1] = {}
        if not id2 in entries: entries[id2] = {}
        value = row[CDTF.VALUE]
        if value == "NA":continue
        entries[id1][id2] = value
        entries[id2][id1] = value
        if not id1 in _s:
            _s.add(id1)
            snps.append(id1)

    snps = list(snps)
    rows = []
    for snp_i in snps:
        row = []
        rows.append(row)
        for snp_j in snps:
            row.append(entries[snp_i][snp_j])
    covariance_matrix = numpy.matrix(rows)
    return snps, covariance_matrix

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



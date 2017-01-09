import pandas
import numpy
import Exceptions

def load_matrix_manager(path):
    d = pandas.read_table(path, sep="\s+")
    m = MatrixManager(d)
    return m


class MatrixManager(object):
    def __init__(self, d):
        _validate(d)
        self.data = _build_data(d)

    def get(self, gene, snps=None):
        return _get(self.data, gene, snps)

    def n_snps(self,gene):
        snps = self.data[gene]
        snps = _non_na(snps)
        snps = {x[CDTF.RSID1] for x in snps}
        return len(snps)

class CDTF(object):
    GENE=0
    RSID1=1
    RSID2=2
    VALUE=3

    K_GENE = "GENE"
    K_RSID1 = "RSID1"
    K_RSID2 = "RSID2"
    K_VALUE = "VALUE"

def _validate(d):
    processed_genes = set()
    last_gene = None
    genes = d[CDTF.K_GENE]
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

def _build_data(d):
    d = d.fillna("NA")
    d.GENE = pandas.Categorical(d.GENE, d.GENE.drop_duplicates())  # speed things up!
    d = zip(d[CDTF.K_GENE].values, d[CDTF.K_RSID1].values, d[CDTF.K_RSID2].values, d[CDTF.K_VALUE].values)
    r = {}
    for t in d:
        gene = t[0]
        if not gene in r:
            r[gene] = []
        r[gene].append(t)
    return r

def _get(d, gene, snps_whitelist=None):
    d = d[gene]

    if snps_whitelist is not None:
        g, r1, r2, v = zip(*d)
        snps = set(r1)
        snps_whitelist = set(snps_whitelist)
        extra = {x for x in snps_whitelist if not x in snps}
        if len(extra):
            msg = "SNPs in whitelist not in matrix:%s"%(extra)
            raise Exceptions.InvalidArguments(msg)
        d = [x for x in d if x[CDTF.RSID1] in snps_whitelist]

    _s = set()
    snps = []
    entries = {}
    for row in d:
        rsid1 = row[CDTF.RSID1]
        rsid2 = row[CDTF.RSID2]
        if not rsid1 in entries: entries[rsid1] = {}
        if not rsid2 in entries: entries[rsid2] = {}
        value = row[CDTF.VALUE]
        if value == "NA":continue
        entries[rsid1][rsid2] = value
        entries[rsid2][rsid1] = value
        if not rsid1 in _s:
            _s.add(rsid1)
            snps.append(rsid1)

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

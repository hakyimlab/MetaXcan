import pandas
import numpy
import Exceptions

GENE="GENE"
RSID1="RSID1"
RSID2="RSID2"
VALUE="VALUE"

class MatrixManager(object):
    def __init__(self, path):
        d = pandas.read_table(path, sep="\s+")
        _validate(d, path)
        self.data = d

    def get(self, gene, snps=None):
        return _get(self.data, gene, snps)


def _get(d, gene, snps_whitelist=None):
    d = d[d[GENE] == gene]

    if snps_whitelist:
        snps = set(d[RSID1])
        snps_whitelist = pandas.Series(snps_whitelist)
        extra = ~ snps_whitelist.isin(snps)
        if numpy.any(extra):
            extra_snps = set(snps_whitelist[extra])
            msg = "SNPs in whitelist not in matrix:%s"%(extra_snps)
            raise Exceptions.InvalidArguments(msg)
        d = d[d[RSID1].isin(snps_whitelist)]

    d = d.dropna()
    snps = d[RSID1].drop_duplicates()
    snps = list(snps)

    entries = {}
    for row in d.iterrows():
        rsid1 = row[1][RSID1]
        rsid2 = row[1][RSID2]
        if not rsid1 in entries: entries[rsid1] = {}
        if not rsid2 in entries: entries[rsid2] = {}
        value = row[1][VALUE]
        entries[rsid1][rsid2] = value
        entries[rsid2][rsid1] = value

    rows = []
    for snp_i in snps:
        row = []
        rows.append(row)
        for snp_j in snps:
            row.append(entries[snp_i][snp_j])
    covariance_matrix = numpy.matrix(rows)
    return snps, covariance_matrix



def _validate(d, path):
    processed_genes = set()
    last_gene = None
    genes = d[GENE]
    for g in genes:
        if g != last_gene:
            if g in processed_genes:
                msg = "Snp Matrix Entries for genes must be contiguous but %s was found in two different places in the file %s." % (g, path)
                raise Exceptions.InvalidInputFormat(msg)
            processed_genes.add(g)
            last_gene = g

    if numpy.any(d.duplicated()):
        raise Exceptions.InvalidInputFormat("Duplicated SNP entries found at %s"%path)
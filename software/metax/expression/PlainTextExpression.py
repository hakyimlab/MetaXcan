import re
import logging
import os
import gzip
import pandas

########################################################################################################################
class ExpressionManager(object):
    def __init__(self, gene_map, file_map):
        self.gene_map = gene_map
        self.file_map = file_map
        self.d = {name:pandas.read_table(path) for name,path in self.file_map.iteritems()}

    def expression_for_gene(self, gene):
        m_ = self.gene_map[gene]
        d = {model:self.d[model][gene].values for model in sorted(m_.keys())}
        return d

    def get_genes(self):
        return self.gene_map.keys()

class ExpressionManagerMemoryEfficient(object):
    def __init__(self, gene_map, file_map):
        self.gene_map = gene_map
        self.file_map = file_map

    def expression_for_gene(self, gene):
        m_ = self.gene_map[gene]
        k = {model:pandas.read_table(self.file_map[model], usecols=[gene])[gene].values for model in sorted(m_.keys())}
        return k

    def get_genes(self):
        return self.gene_map.keys()

_exregex = re.compile("TW_(.*)_0.5.expr.txt")
def _structure(folder, pattern=None):
    logging.info("Acquiring expression files")
    files = os.listdir(folder)
    gene_map = {}
    file_map = {}

    _regex = _exregex if not pattern else re.compile(pattern)

    for file in files:
        if _regex.search(file):
            name = _regex.match(file).group(1)
            name = name.replace("-","_")
        else:
            continue
        path = os.path.join(folder, file)
        _o = gzip.open if ".gz" in file else open
        with _o(path) as f:
            comps = f.readline().strip().split()

            for i,gene in enumerate(comps):
                if not gene in gene_map:
                    gene_map[gene] = {}
                gene_map[gene][name] = i

            file_map[name] = path

    return gene_map, file_map

########################################################################################################################

class Expression(object):
    def __init__(self, path):
        self.d = pandas.read_table(path)

    def expression_for_gene(self, gene):
        k = self.d[gene].values
        return k

    def get_genes(self):
        return self.d.columns.values
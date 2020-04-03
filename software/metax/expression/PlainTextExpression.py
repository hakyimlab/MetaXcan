import re
import logging
import os
import io
import gzip
import pandas

from . import Expression
from ..misc import Math

########################################################################################################################
class ExpressionManager(Expression.ExpressionManager):
    def __init__(self, folder, pattern=None, standardise=False):
        gene_map, file_map = _structure(folder, pattern)
        self.gene_map = gene_map
        self.file_map = file_map
        self.d = None
        self.standardise = standardise

    def expression_for_gene(self, gene):
        m_ = self.gene_map[gene]
        d = {model:self.d[model][gene].values for model in sorted(m_.keys())}

        if self.standardise:
            k_ = {}
            for key, value in d.items():
                t_ = Math.standardize(value)

                if t_ is not None:
                    k_[key] = t_
            d = k_
        return d

    def get_genes(self):
        return list(self.gene_map.keys())

    def enter(self):
        d = {}
        for name in sorted(self.file_map.keys()):
            path = self.file_map[name]
            logging.log(9, "Loading %s", path)
            d_ = pandas.read_table(path)
            if "FID" in d_:
                d_ = d_.drop(["FID", "IID"], axis=1)
            d[name] = d_
        self.d = d

    def exit(self):
        pass

class ExpressionManagerMemoryEfficient(Expression.ExpressionManager):
    def __init__(self, folder, pattern, standardise = False):
        gene_map, file_map = _structure(folder, pattern)
        self.gene_map = gene_map
        self.file_map = file_map
        self.standardise = standardise

    def expression_for_gene(self, gene):
        m_ = self.gene_map[gene]
        k = {model:pandas.read_table(self.file_map[model], usecols=[gene])[gene].values for model in sorted(m_.keys())}
        if self.standardise:
            k_ = {}
            for key, value in k.items():
                t_ = Math.standardize(value)
                if t_ is not None:
                    k_[key] = t_
            k = k_
        return k

    def get_genes(self):
        return list(self.gene_map.keys())

_exregex = re.compile("TW_(.*)_0.5.expr.txt")
def _structure(folder, pattern=None):
    logging.info("Acquiring expression files")
    files = os.listdir(folder)
    gene_map = {}
    file_map = {}

    _regex = re.compile(pattern) if pattern else None

    for file in files:
        if _regex:
            if _regex.search(file):
                name = _regex.match(file).group(1)
            else:
                continue
        else:
            name = file
        name = name.replace("-", "_")
        path = os.path.join(folder, file)

        def _ogz(p):
            return io.TextIOWrapper(gzip.open(p, "r"), newline="")
        _o = _ogz if ".gz" in path else open
        with _o(path) as f:
            comps = f.readline().strip().split()

            for i,gene in enumerate(comps):
                if gene == "FID" or gene == "IID":
                    continue
                if not gene in gene_map:
                    gene_map[gene] = {}
                gene_map[gene][name] = i

            file_map[name] = path

    return gene_map, file_map

########################################################################################################################

class Expression(Expression.Expression):
    def __init__(self, path):
        self.path = path
        self.d = None

    def expression_for_gene(self, gene):
        k = self.d[gene].values
        return k

    def get_genes(self):
        return self.d.columns.values

    def enter(self):
        self.d = pandas.read_table(self.path)
        if "FID" in self.d:
            self.d = self.d.drop(["FID", "IID"], axis = 1)

    def exit(self):
        pass
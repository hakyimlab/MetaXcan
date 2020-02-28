import os
import re
import numpy
import logging
import h5py
try:
    import h5py_cache
except:
    logging.info("Couldn't import h5py_cache. Anyway, this dependency should be removed. It has been folded into h5py")

from . import Expression as _Expression

from ..misc import Math

########################################################################################################################

class ExpressionManager(_Expression.ExpressionManager):
    def __init__(self, folder, pattern=None, code_999=False, standardise=False):
        self.gene_map = None
        self.h5 = None
        self.folder = folder
        self.pattern = pattern
        self.code_999 = code_999
        self.last_gene = None
        self.k = None
        self.standardise = standardise

    def expression_for_gene(self, gene):
        if gene == self.last_gene:
            return self.k

        m_ = self.gene_map[gene]
        k = {model:self.h5[model]['pred_expr'][m_[model]] for model in sorted(m_.keys())}
        if self.code_999:
            logging.log(8, "Looking for code -999")
            k = _code_999(k)

        if self.standardise:
            logging.log(8, "Standardizing")
            k_ = {}
            for key, value in k.items():
                t_ = Math.standardize(value)
                if t_ is not None:
                    k_[key] = t_
            k = k_


        self.k = k
        self.last_gene = gene

        return k

    def get_genes(self):
        return list(self.gene_map.keys())

    def enter(self):
        gene_map, h5 = _structure(self.folder, self.pattern)
        self.gene_map = gene_map
        self.h5 = h5
        return self

    def exit(self):
        _close(self.h5)

def _structure(folder, pattern=None):
    logging.info("Acquiring HDF5 expression files")
    files = os.listdir(folder)
    h5 = {}
    gene_map = {}

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
        file = h5py.File(path, 'r')

        genes = [g for g in file['genes']]
        for i,gene in enumerate(genes):
            if not gene in gene_map:
                gene_map[gene] = {}
            gene_map[gene][name] = i

        h5[name] = file

    return gene_map, h5

def _close(h5):
    for name, h_ in h5.items():
        h_.close()

def _code_999(k):
    logged_ = []

    def c_(x):
        close = numpy.isclose(x, -999.0, atol=1e-3, rtol=0)
        if numpy.sum(close):
            if not logged_:
                logged_.append(True)
                logging.warning("Expression Manager: Encountered value of -999")
            x[close] = numpy.nan
        return x

    k = {k:c_(v) for k, v in k.items()}
    return  k

########################################################################################################################

class Expression(_Expression.Expression):
    def __init__(self, path, code_999=False):
        self.path = path
        self.genes = None
        self.h5 = None
        self.gene_idx = None
        self.code_999 = code_999

    def expression_for_gene(self, gene):
        k = self.h5[self.gene_idx[gene]]
        if self.code_999:
            k = _code_999_b(k)
        return k

    def get_genes(self):
        return self.genes

    def enter(self):
        genes, h5 = _structure_file(self.path)
        self.h5 = h5
        self.genes = genes
        self.gene_idx = {gene: id for id, gene in enumerate(genes)}

    def exit(self):
        _close_file(self.h5)

def _code_999_b(k):
    logged_ = []

    def c_(x):
        if numpy.isclose(x, -999, atol=1e-3, rtol=0):
            if not logged_:
                logged_.append(True)
                logging.warning("Expression: Encountered value of -999")
            return numpy.nan
        return x

    k = numpy.array(list(map(c_, k)))
    return  k

def _structure_file(file_path):
    logging.info("Acquiring HDF5 expression cache")
    file = h5py.File(file_path, 'r')
    genes = [g for g in file['genes']]
    h5 = file['pred_expr']
    return genes, h5

def _close_file(h5):
    h5.file.close()

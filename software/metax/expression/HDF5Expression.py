import os
import re
import numpy
import logging
import h5py
import h5py_cache

import Expression as _Expression

########################################################################################################################

class ExpressionManager(_Expression.ExpressionManager):
    def __init__(self, folder, pattern=None, code_999=False):
        self.gene_map = None
        self.h5 = None
        self.folder = folder
        self.pattern = pattern
        self.code_999 = code_999

    def expression_for_gene(self, gene):
        m_ = self.gene_map[gene]
        k = {model:self.h5[model]['pred_expr'][m_[model]] for model in sorted(m_.keys())}
        if self.code_999:
            k = _code_999(k)
        return k

    def get_genes(self):
        return self.gene_map.keys()

    def enter(self):
        gene_map, h5 = _structure(self.folder, self.pattern)
        self.gene_map = gene_map
        self.h5 = h5

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
    for name, h_ in h5.iteritems():
        h_.close()

def _code_999(k):
    logged_ = []

    def c_(x):
        if numpy.isclose(x, -999, atol=1e-3, rtol=0):
            if not logged_:
                logged_.append(True)
                logging.warning("Expression Manager: Encountered value of -999")
            return numpy.nan
        return x

    k = {k: numpy.array(map(c_, v)) for k, v in k.iteritems()}
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

    k = numpy.array(map(c_, k))
    return  k

def _structure_file(file_path):
    logging.info("Acquiring HDF5 expression cache")
    HDF5_CACHE = int(60 * (1024 ** 2))
    file =  h5py_cache.File(file_path, 'r', chunk_cache_mem_size=HDF5_CACHE, w0=1.0, dtype='float32')
    genes = [g for g in file['genes']]
    h5 = file['pred_expr']
    return genes, h5

def _close_file(h5):
    h5.file.close()

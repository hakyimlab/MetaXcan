import os
import re
import h5py_cache

class ExpressionManager(object):
    def __init__(self, gene_map, h5):
        self.gene_map = gene_map
        self.h5 = h5
        self.own = False

_exregex = re.compile("pred_TW_(.*)_0.5_hrc_hapmap.h5")
def _structure(folder):
    files = os.listdir(folder)
    h5 = {}
    gene_map = {}

    HDF5_CACHE = int(500 * (1024 ** 2))

    for file in files:
        if _exregex.search(file):
            name = _exregex.match(file).group(1)
        else:
            continue
        path = os.path.join(folder, file)
        file = h5py_cache.File(path, 'r', chunk_cache_mem_size=HDF5_CACHE)

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
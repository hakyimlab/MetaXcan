import logging
import os

import numpy
import pandas

from .. import NamingConventions
from .. import Utilities


class MetaXcanResultsManager(object):
    def __init__(self, data):
        data, genes, models = _build_data(data)
        self.data = data
        self.genes = genes
        self.models = models

    def results_for_gene(self, gene):
        if not gene in self.data: return None, None
        results = self.data[gene]
        labels = sorted(results.keys())
        values = [results[k] for k in labels]
        return values, labels

    def get_genes(self):
        return self.genes

    def get_model_labels(self):
        return self.models

def build_manager(folder, filters=[".*csv"]):
    files = Utilities.target_files(folder, file_filters=filters)
    results = _load_results(files)
    manager = MetaXcanResultsManager(results)
    return manager

def _load_results(files):
    results = {}
    for file in sorted(files):
        logging.log(9, "Loading metaxcan %s", file)
        e = pandas.read_csv(file)
        root, name = os.path.split(file)
        pheno, model, model_tag, model_type  = NamingConventions.parse_file_name(name)
        e["pheno"] = pheno
        e["tissue"] = model
        results[model] = e
    return results

def _build_data(data):
    logging.log(9,"Building data")
    genes = set()
    models = set()
    result = {}
    for k, df in data.iteritems():
        logging.log(9, "Processing %s", k)
        for i, row in df.iterrows():
            if not numpy.isfinite(row.zscore):
                continue
            if not row.gene in result:
                result[row.gene]={}
            gene_entries = result[row.gene]
            gene_entries[row.tissue] = row.zscore
            genes.add(row.gene)
            models.add(row.tissue)
    return result, genes, models

def _get_columns(data):
    logging.info("Getting columns")
    columns = {x for k in data.keys() for x in data[k].columns}
    return columns
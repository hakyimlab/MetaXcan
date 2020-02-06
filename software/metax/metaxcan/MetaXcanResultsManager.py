import logging
import os
import re

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

def build_manager(folder, filters=[".*csv"], file_name_pattern=None):
    files = Utilities.target_files(folder, file_filters=filters)
    results = _load_results(files, file_name_pattern)
    manager = MetaXcanResultsManager(results)
    return manager

def _parse_name(name, file_name_pattern):
    if not file_name_pattern:
        pheno, model, model_tag, model_type = NamingConventions.parse_file_name(name)
        return pheno, model

    r = re.compile(file_name_pattern)
    g = r.match(name).groups()
    if len(g) > 1:
        pheno, model = g[0], g[1]
    else:
        pheno, model = "", g[0]
    return pheno, model

def _load_results(files, file_name_pattern):
    results = {}
    for file in sorted(files):
        logging.log(9, "Loading metaxcan %s", file)
        e = pandas.read_csv(file)
        root, name = os.path.split(file)
        pheno, model = _parse_name(name, file_name_pattern)
        e["pheno"] = pheno
        e["tissue"] = model
        results[model] = e
    return results

def _build_data(data):
    logging.log(9,"Building data")
    genes = set()
    models = set()
    result = {}
    for k, df in data.items():
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
    columns = {x for k in list(data.keys()) for x in data[k].columns}
    return columns
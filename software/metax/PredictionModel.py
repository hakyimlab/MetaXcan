import os
import sqlite3
import pandas
import logging

import Exceptions

class WDBQF(object):
    "Weight DB weight Query Format"
    RSID=0
    GENE=1
    WEIGHT=2
    REF_ALLELE=3
    EFF_ALLELE=4

    K_RSID="rsid"
    K_GENE="gene"
    K_WEIGHT="weight"
    K_EFFECT_ALLELE="effect_allele"
    K_NON_EFFECT_ALLELE="non_effect_allele"

    ORDER = [(K_RSID, RSID), (K_GENE,GENE), (K_WEIGHT,WEIGHT), (K_EFFECT_ALLELE, EFF_ALLELE),(K_NON_EFFECT_ALLELE, REF_ALLELE)]

class WDBEQF(object):
    "Weight DB extra table Query Format"
    GENE=0
    GENE_NAME=1
    N_SNP_IN_MODEL=2
    PRED_PERF_R2=3
    PRED_PERF_PVAL=4
    PRED_PERF_QVAL=5

    K_GENE="gene"
    K_GENE_NAME="gene_name"
    K_N_SNP_IN_MODEL="n_snps_in_model"
    K_PRED_PERF_R2="pred_perf_r2"
    K_PRED_PERF_PVAL="pred_perf_pval"
    K_PRED_PERF_QVAL="pred_perf_qval"

    ORDER =[(K_GENE,GENE), (K_GENE_NAME,GENE_NAME), (K_N_SNP_IN_MODEL, N_SNP_IN_MODEL), (K_PRED_PERF_R2,PRED_PERF_R2), (K_PRED_PERF_PVAL, PRED_PERF_PVAL), (K_PRED_PERF_QVAL, PRED_PERF_QVAL)]

def snps_in_db(path):
    t = ModelDB(path)
    w = t.load_weights()
    s = set(w[WDBQF.RSID])
    return s

class ModelDB(object):
    def __init__(self, file_name , create_if_absent=False):
        self.connection = None
        self.cursor = None
        self.file_name = file_name
        self.create_if_absent = create_if_absent

    def __del__(self):
        self.closeDB()

    def openDBIfNecessary(self):
        if not self.connection:
            if not self.create_if_absent and not os.path.exists(self.file_name):
                raise Exceptions.BadFilename(self.file_name)
            self.connection = sqlite3.connect(self.file_name)
            self.cursor = self.connection.cursor()

    def closeDB(self):
        if self.connection:
            self.connection.close()
            self.connection = None
            self.cursor = None

    def load_weights(self, gene_key=None):
        self.openDBIfNecessary()

        params = []
        query, params = query_helper("SELECT rsid, gene, weight, ref_allele, eff_allele FROM weights", gene_key)

        try:
            results = self.cursor.execute(query, params)
        except sqlite3.OperationalError as e:
            raise Exceptions.ReportableException("Could not read input tissue database. Please try updating the tissue model files.")
        except Exception as e:
            raise e

        weights = zip(*results)
        return  weights

    def load_extra(self, gene_key=None):
        self.openDBIfNecessary()

        query,params = query_helper("SELECT gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval` FROM extra", gene_key)

        try:
            results = self.cursor.execute(query, params)
        except sqlite3.OperationalError as e:
            print str(e)
            raise Exceptions.ReportableException("Could not read input tissue database. Please try updating the tissue model files.")
        except Exception as e:
            raise e

        extra = zip(*results)
        return  extra

def query_helper(query, gene_key=None):
    params = []
    if gene_key:
        query += " WHERE gene = ?"
        params.append(gene_key)
    query += ";"
    return query, tuple(params)

class Model(object):
    def __init__(self, weights, extra):
        self.weights = weights
        self.extra = extra

    def snps(self):
        snps = self.weights.rsid.values
        return set(snps)

def dataframe_from_weight_data(w):
    weights = pandas.DataFrame({key: w[order] for key, order in WDBQF.ORDER})
    weights = weights[[key for key,order in WDBQF.ORDER]]
    return weights

def dataframe_from_extra_data(e):
    extra = pandas.DataFrame({key: e[order] for key, order in WDBEQF.ORDER})
    extra = extra[[key for key, order in WDBEQF.ORDER]]
    return extra

def load_model(path):
    db = ModelDB(path)
    weights, extra = db.load_weights(), db.load_extra()
    weights = dataframe_from_weight_data(weights)
    extra = dataframe_from_extra_data(extra)
    model = Model(weights, extra)
    return model

class ModelManager(object):
    def __init__(self, models):
        self.models = models #models is a dictionary of dictionaries, models[gene][model] are the weights for a given gene and model

    def get_genes(self):
        return set(self.models.index.get_level_values(0))

    def get_snps(self, gene = None):
        w = self.models if not gene else self.models.loc[gene]
        i = 2 if not gene else 1
        snps = set(w.index.get_level_values(i))
        return snps

    def get_weights(self, gene):
        return self.models.loc[gene]

def _split_models(models):
    _m = {}
    for k,m in models.iteritems():
        for gene in m.extra.gene:
            w = m.weights
            w = w[w.gene == gene]
            if not gene in _m: _m[gene] = {}
            _m[gene][k] = w
    return _m


def load_model_manager(path):
    def _extract_model_name(path):
        p = os.path.split(path)[1]
        p = p.split("_0.5.db")[0]
        return p

    def _get_models(paths):
        logging.log(9, "preloading models")
        _m = {_extract_model_name(x): load_model(x) for x in paths}
        for k,m in _m.iteritems():
            logging.log(9, "processing %s", k)
            w = m.weights
            w["model"] = k
        _m = [x.weights for x in _m.values()]
        models = pandas.concat(_m)
        logging.log(9, "indexing")
        models = models.set_index(["gene", "model", "rsid"])
        models = models.sort_index()
        return models

    paths = [os.path.join(path,x) for x in os.listdir(path) if ".db" in x]
    models = _get_models(paths)
    model_manager = ModelManager(models)
    return model_manager

import logging
import os
import sqlite3
import re
import pandas

from . import Exceptions
from . import NamingConventions


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
    def __init__(self, file_name, create_if_absent=False, snp_key=None):
        self.connection = None
        self.cursor = None
        self.file_name = file_name
        self.create_if_absent = create_if_absent
        self.snp_key = snp_key

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
        _q = "SELECT {}, gene, weight, ref_allele, eff_allele FROM weights".format(self.snp_key if self.snp_key else "rsid")
        query, params = query_helper(_q, gene_key)

        try:
            results = self.cursor.execute(query, params)
        except sqlite3.OperationalError as e:
            raise Exceptions.ReportableException("Could not read input tissue database. Please try updating the tissue model files.")
        except Exception as e:
            raise e

        weights = list(zip(*results))
        return  weights

    def load_extra(self, gene_key=None):
        self.openDBIfNecessary()

        query,params = query_helper("SELECT gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval` FROM extra order by gene", gene_key)

        try:
            results = self.cursor.execute(query, params)
        except sqlite3.OperationalError as e:
            logging.info(str(e))
            raise Exceptions.ReportableException("Could not read input tissue database. Please try updating the tissue model files.")
        except Exception as e:
            raise e

        extra = list(zip(*results))
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
    if len(w) == 0:
        w = [[] for key, order in WDBQF.ORDER]
    weights = pandas.DataFrame({key: w[order] for key, order in WDBQF.ORDER})
    weights = weights[[key for key,order in WDBQF.ORDER]]
    return weights

def dataframe_from_extra_data(e):
    if len(e) == 0:
        e = [[] for key, order in WDBEQF.ORDER]
    extra = pandas.DataFrame({key: e[order] for key, order in WDBEQF.ORDER})
    extra = extra[[key for key, order in WDBEQF.ORDER]]
    return extra

def load_model(path, snp_key=None):
    db = ModelDB(path, snp_key=snp_key)
    weights, extra = db.load_weights(), db.load_extra()

    #UGLY patch. UGLY. UGLY.
    if snp_key:
        s = weights[WDBQF.RSID]
        if "_b38" in s[0] and not "chr" in s[0]:
            logging.info("Patching variant names with -chr- prefix")
            s = ["chr"+x for x in s]
            weights[WDBQF.RSID] = s

    weights = dataframe_from_weight_data(weights)
    extra = dataframe_from_extra_data(extra)

    model = Model(weights, extra)

    return model

def load_genes(folder, name_filter=None):
    model_paths = _model_paths(folder, name_filter)
    models = pandas.DataFrame()

    for path in model_paths:
        m = load_model(path)
        e = m.extra
        e = e[["gene", "gene_name"]]
        models = pandas.concat([models, e])

    models = models.drop_duplicates()
    return models

###############################################################################
class ModelManagerBase(object):
    def get_genes(self): raise Exceptions.NotImplemented("ModelManager: get_genes")
    def get_implicated_genes(self, snps): raise Exceptions.NotImplemented("ModelManager: get_implicated_genes")
    def get_rsids(self, gene = None):raise Exceptions.NotImplemented("ModelManager: get_rsids")
    def get_model_labels(self, gene = None):raise Exceptions.NotImplemented("ModelManager: get_model_labels")
    def get_models(self, gene): raise Exceptions.NotImplemented("ModelManager: get_models")

###############################################################################
class ModelManager(ModelManagerBase):
    def __init__(self, models):
        models = _prepare_models(models)
        self.models = models
        self.snp_keys = _get_snp_key(models)

    def get_genes(self):
        return set(self.models.index.get_level_values(0))

    def get_implicated_genes(self, snps):
        return  _get_implicated(self.snp_keys, snps)

    def get_rsids(self, gene = None):
        w = self.models if not gene else self.models.loc[gene]
        i = 2 if not gene else 1
        snps = set(w.index.get_level_values(i))
        return snps

    def get_model_labels(self, gene = None):
        w = self.models if not gene else self.models.loc[gene]
        i = 1 if not gene else 0
        labels = set(w.index.get_level_values(i))
        return labels

    def get_models(self, gene):
        return self.models.loc[gene]

def _get_implicated(snp_keys, snps):
    genes = set()
    for snp in snps:
        if snp in snp_keys:
            genes.update(snp_keys[snp])
    return genes

def _prepare_models(models):
    logging.log(9, "preparing models (indexing)")
    models = models.set_index(["gene", "model", "rsid"])
    models = models.sort_index()
    return models

def _get_snp_key(models):
    logging.log(9, "preparing snp keys")
    keys = {}
    for k in models.itertuples():
        rsid = k.Index[2]
        gene = k.Index[0]
        if not rsid in keys: keys[rsid] = set()
        keys[rsid].add(gene)
    return keys

def _model_paths(path, name_filter=None):
    f = [re.compile(x) for x in name_filter] if name_filter else [re.compile(".*db$")]
    paths = [os.path.join(path, x) for x in os.listdir(path) if True in [f_.search(x) is not None for f_ in f]]
    return paths

###############################################################################
class _ModelManager(ModelManagerBase):
    """Version that performs certain operations faster, but returns data in different format!"""
    def __init__(self, models):
        models, rsids, snp_key = _prepare_models_2(models)
        self.models = models
        self.rsids = rsids
        self.snp_key = snp_key

    def get_genes(self):
        return set(self.models.keys())

    def get_implicated_genes(self, snps):
        return  _get_implicated(self.snp_keys, snps)

    def get_rsids(self, gene = None):
        if not gene: return self.rsids
        if not gene in self.models: return None

        g = self.models[gene]
        rsids = set()
        for tissue, weights in g.items():
            rsids.update(list(weights.keys()))
        return rsids

    def get_model_labels(self, gene = None):
        if not gene:
            labels = set()
            for gene,tissues in self.models.items():
                labels.update(list(tissues.keys()))
            return labels

        if not gene in self.models: return None
        return set(self.models[gene].keys())

    def get_models(self, gene):
        return self.models[gene]

def _prepare_models_2(models):
    logging.log(9, "preparing models (dictionary layout)")
    rsids = set()
    rsid_to_genes = {}

    r = {}
    for t in models.itertuples():
        if not t.gene in r: r[t.gene] = {}
        g = r[t.gene]
        if not t.model in g: g[t.model] = {}
        m = g[t.model]
        m[t.rsid] = t.weight
        rsids.add(t.rsid)

        if not t.rsid in rsid_to_genes: rsid_to_genes[t.rsid] = set()
        rsid_to_genes[t.rsid].add(t.gene)
    return r, rsids, rsid_to_genes

###############################################################################
def load_model_manager(path, trim_ensemble_version=False, Klass=ModelManager, name_pattern=None, name_filter=None, model_db_snp_key=None):

    def _get_models(paths, trim_ensemble_version=False):
        logging.log(9, "preloading models")
        _m = {NamingConventions.extract_model_name(x, name_pattern): load_model(x, model_db_snp_key) for x in paths}
        keys = sorted(_m.keys())
        for i,k in enumerate(keys):
            logging.log(9, "processing %s", k)
            m = _m[k]
            w = m.weights
            w["model"] = k
        _m = [x.weights for x in list(_m.values())]
        models = pandas.concat(_m)
        if trim_ensemble_version:
            k = models.gene.str.split(".").str.get(0)
            if len(set(k)) != len(set(models.gene)):
                raise Exceptions.ReportableException("genes cannot lose the ensemble version id")
            models.gene = k
        return models

    paths = _model_paths(path, name_filter)
    models = _get_models(paths, trim_ensemble_version)
    model_manager = Klass(models)
    return model_manager

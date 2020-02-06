__author__ = 'heroico'

import sqlite3
import os
from collections import OrderedDict
from . import Exceptions

class GeneEntry:
    def __init__(self, gene, gene_name, n_snps, R2, pval,qval):
        self.gene = gene
        self.gene_name = gene_name
        self.n_snps = n_snps
        self.pred_perf_R2 = R2
        self.pred_perf_pval = pval
        self.pred_perf_qval = qval

class WeightDBEntry:
    def __init__(self, rsid=None, gene=None, weight=None, ref_allele=None, eff_allele=None, pval=None, N=None, cis=None):
        """Warning: many db's have empty 'N', 'cis' and 'pval'"""
        self.rsid = rsid
        self.gene = gene
        self.weight = weight
        self.ref_allele = ref_allele
        self.eff_allele = eff_allele
        self.pval = pval
        self.N = N
        self.cis = cis

class WDBQF(object):
    "Weight DB weight Query Format"
    RSID=0
    GENE=1
    WEIGHT=2
    REF_ALLELE=3
    EFF_ALLELE=4

class WDBEQF(object):
    "Weight DB extra table Query Format"
    GENE=0
    GENE_NAME=1
    N_SNP_IN_MODEL=2
    PRED_PERF_R2=3
    PRED_PERF_PVAL=4
    PRED_PERF_QVAL=5

class WeightDB(object):
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
                raise RuntimeError("Weight file doesn't exist")
            self.connection = sqlite3.connect(self.file_name)
            self.cursor = self.connection.cursor()

    def closeDB(self):
        if self.connection:
            self.connection.close()
            self.connection = None
            self.cursor = None

    def weightEntriesFromResults(self, results, extra, result_callback=None):
        weights = []
        for result in results:
            weight = WeightDBEntry(result[WDBQF.RSID],
                                   result[WDBQF.GENE],
                                   result[WDBQF.WEIGHT],
                                   result[WDBQF.REF_ALLELE],
                                   result[WDBQF.EFF_ALLELE])
            weights.append(weight)
            if result_callback:
                result_callback(weight, extra)
        return weights

    def loadFromDB(self, callback=None, gene_key=None):
        self.openDBIfNecessary()
        extra = self.loadExtraColumnData()
        extra = {e.gene:e for e in extra}

        if gene_key is None:
            results = self.cursor.execute("SELECT rsid, gene, weight, ref_allele, eff_allele FROM weights;")
        else:
            results = self.cursor.execute("SELECT rsid, gene, weight, ref_allele, eff_allele FROM weights where gene = ?;", (gene_key))

        weights = self.weightEntriesFromResults(results, extra, callback)
        return  weights

    def loadExtraColumnData(self, gene_key=None):
        self.openDBIfNecessary()
        try:
            if gene_key is None:
                results = self.cursor.execute("SELECT gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval` FROM extra;")
            else:
                results = self.cursor.execute("SELECT gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval` FROM extra WHERE gene = ?;", (gene_key,))
        except sqlite3.OperationalError as e:
            print(str(e))
            raise Exceptions.ReportableException("Could not read input tissue database. Please try updating the tissue model files.")
        except Exception as e:
            raise e

        extra = [GeneEntry(x[WDBEQF.GENE], x[WDBEQF.GENE_NAME], x[WDBEQF.N_SNP_IN_MODEL], x[WDBEQF.PRED_PERF_R2], x[WDBEQF.PRED_PERF_PVAL], x[WDBEQF.PRED_PERF_QVAL]) for x in results]
        return  extra

    def loadGeneNamesFromDB(self):
        self.openDBIfNecessary()
        names = []

        results = self.cursor.execute("SELECT DISTINCT gene FROM weights;")

        for result in results:
            name = result[0]
            names.append(name)

        return names

class WeightDBEntryLogic(object):
    def __init__(self, db_file_name):
        self.weights_by_gene = OrderedDict()#{}
        self.genes_for_an_rsid = OrderedDict()#{}
        self.gene_data_for_gene = OrderedDict()#{}
        self._loadData(db_file_name)

    def anEntryWithRSID(self, rsid):
        entry = None
        if not rsid in self.genes_for_an_rsid:
            return entry

        genes = self.genes_for_an_rsid[rsid]
        gene = genes[0]
        weights = self.weights_by_gene[gene]
        entry = weights[rsid]
        return  entry

    def _loadData(self, db_file_name):
        weights_db = WeightDB(db_file_name)

        class ByNameCallback(object):
            """Helper class to group weights by gene name"""
            def __init__(self, weights_by_gene, genes_for_an_rsid, gene_data_for_gene):
                self.weights_by_gene = weights_by_gene
                self.genes_for_an_rsid = genes_for_an_rsid
                self.gene_data_for_gene = gene_data_for_gene

            def __call__(self, weight, extra):
                if weight.gene in self.weights_by_gene:
                    weights = self.weights_by_gene[weight.gene]
                else:
                    weights = OrderedDict()
                    self.weights_by_gene[weight.gene] = weights
                weights[weight.rsid]= weight

                if not weight.rsid in self.genes_for_an_rsid:
                    self.genes_for_an_rsid[weight.rsid] = []
                genes = self.genes_for_an_rsid[weight.rsid]

                if not weight.gene in genes:
                    genes.append(weight.gene)

                gene_entry = extra[weight.gene]
                self.gene_data_for_gene[weight.gene] = gene_entry

        callback = ByNameCallback(self.weights_by_gene, self.genes_for_an_rsid, self.gene_data_for_gene)
        weights_db.loadFromDB(callback)

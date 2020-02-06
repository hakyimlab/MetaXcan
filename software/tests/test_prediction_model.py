import os
import numpy
import numpy.testing

import pandas
from sqlalchemy import create_engine

import unittest
from  metax import PredictionModel

from . import SampleData

def get_model_weights(path):
    engine = create_engine('sqlite:///'+path)
    return pandas.read_sql_table('weights', engine)

def get_weights_in_models(path):
    f = [x for x in os.listdir(path) if ".weights" in x]
    r = []
    for f_ in f:
        p_ = os.path.join(path,f_)
        weights = pandas.read_table(p_)
        weights["model"] = f_.split(".weights")[0]
        r.append(weights)
    return pandas.concat(r)

def _compare(unit_test, m, w):
    unit_test.assertEqual(set(w.index.values), set(m.index.values))
    numpy.testing.assert_array_equal(w.loc[m.index.values].ref_allele, m.ref_allele)
    numpy.testing.assert_array_equal(w.loc[m.index.values].eff_allele, m.eff_allele)
    numpy.testing.assert_array_almost_equal(w.loc[m.index.values].weight, m.weight)

def _compare_o(unit_test, m, w):
    for t in w.itertuples():
        numpy.testing.assert_almost_equal(m[t.model][t.rsid], t.weight)

def compare_model_manager_models_to_weights(unit_test,model_manager, weights, model, gene):
    w = weights[(weights.model == model) & (weights.gene == gene)][["rsid", "ref_allele", "eff_allele", "weight"]].set_index("rsid")
    m = model_manager.models.loc[gene, model].rename(columns={"effect_allele":"eff_allele", "non_effect_allele":"ref_allele"})
    _compare(unit_test, m, w)

def _assert_optimized_manager(unit_test, model_manager, weights, genes):
    unit_test.assertEqual(model_manager.get_genes(), genes)
    unit_test.assertEqual(model_manager.get_rsids(), set(weights.rsid))
    unit_test.assertEqual(model_manager.get_model_labels(), set(weights.model))
    #
    for gene in set(weights.gene):
        w = weights[weights.gene == gene]
        m = model_manager.get_models(gene)
        _compare_o(unit_test, m, w)

    #
    for gene in set(weights.gene):
        unit_test.assertEqual(model_manager.get_model_labels(gene), set(weights[weights.gene == gene].model))

class TestPredictionModelDB(unittest.TestCase):

    def test_load(self):
        t = PredictionModel.ModelDB("tests/_td/dbs/test_1.db")
        extra = t.load_extra()
        self.assertEqual(len(extra), 6)

        e_e = list(zip(*(SampleData.sample_extra_2())))
        self.assertEqual(extra[PredictionModel.WDBEQF.GENE], e_e[PredictionModel.WDBEQF.GENE])
        self.assertEqual(extra[PredictionModel.WDBEQF.GENE_NAME], e_e[PredictionModel.WDBEQF.GENE_NAME])
        self.assertEqual(extra[PredictionModel.WDBEQF.N_SNP_IN_MODEL], e_e[PredictionModel.WDBEQF.N_SNP_IN_MODEL])
        self.assertEqual(extra[PredictionModel.WDBEQF.PRED_PERF_R2], e_e[PredictionModel.WDBEQF.PRED_PERF_R2])
        self.assertEqual(extra[PredictionModel.WDBEQF.PRED_PERF_PVAL], e_e[PredictionModel.WDBEQF.PRED_PERF_PVAL])
        self.assertEqual(extra[PredictionModel.WDBEQF.PRED_PERF_QVAL], e_e[PredictionModel.WDBEQF.PRED_PERF_QVAL])

        weights = t.load_weights()
        self.assertEqual(len(weights), 5)

        e_w = list(zip(*(SampleData.sample_weights_2())))
        self.assertEqual(weights[PredictionModel.WDBQF.RSID], e_w[PredictionModel.WDBQF.RSID])
        self.assertEqual(weights[PredictionModel.WDBQF.GENE], e_w[PredictionModel.WDBQF.GENE])
        self.assertEqual(weights[PredictionModel.WDBQF.WEIGHT], e_w[PredictionModel.WDBQF.WEIGHT])
        self.assertEqual(weights[PredictionModel.WDBQF.REF_ALLELE], e_w[PredictionModel.WDBQF.REF_ALLELE])
        self.assertEqual(weights[PredictionModel.WDBQF.EFF_ALLELE], e_w[PredictionModel.WDBQF.EFF_ALLELE])

    def test_snps_in_db(self):
        expected = {"rs245915", "rs245913", "rs245909", "rs245906", "rs10486599", "rs144012121", "rs117887801", "rs542000", "rs544632",
                    "rs498475", "rs849327", "rs849336", "rs849335", "rs1513272", "rs849135", "rs849134", "rs860262", "rs849133", "rs1635852",
                    "rs864745", "rs112751321", "rs144273091", "rs117462481", "rs149305679", "rs643036", "rs1937888", "rs17155745", "rs62626328"}
        actual = PredictionModel.snps_in_db("tests/_td/dbs/test_2.db")
        self.assertEqual(actual, expected)

    def test_load_model(self):
        snp_model = PredictionModel.load_model("tests/_td/dbs/test_1.db")

        e_e = SampleData.dataframe_from_extra(SampleData.sample_extra_2())
        numpy.testing.assert_array_equal(snp_model.extra[PredictionModel.WDBEQF.K_GENE], e_e[PredictionModel.WDBEQF.K_GENE])
        numpy.testing.assert_array_equal(snp_model.extra[PredictionModel.WDBEQF.K_GENE_NAME], e_e[PredictionModel.WDBEQF.K_GENE_NAME])
        numpy.testing.assert_array_equal(snp_model.extra[PredictionModel.WDBEQF.K_N_SNP_IN_MODEL], e_e[PredictionModel.WDBEQF.K_N_SNP_IN_MODEL])
        numpy.testing.assert_array_equal(snp_model.extra[PredictionModel.WDBEQF.K_PRED_PERF_R2], e_e[PredictionModel.WDBEQF.K_PRED_PERF_R2])
        numpy.testing.assert_array_equal(snp_model.extra[PredictionModel.WDBEQF.K_PRED_PERF_PVAL], e_e[PredictionModel.WDBEQF.K_PRED_PERF_PVAL])
        numpy.testing.assert_array_equal(snp_model.extra[PredictionModel.WDBEQF.K_PRED_PERF_QVAL], e_e[PredictionModel.WDBEQF.K_PRED_PERF_QVAL])

        e_w = SampleData.dataframe_from_weights(SampleData.sample_weights_2())
        numpy.testing.assert_array_equal(snp_model.weights[PredictionModel.WDBQF.K_RSID], e_w[PredictionModel.WDBQF.K_RSID])
        numpy.testing.assert_array_equal(snp_model.weights[PredictionModel.WDBQF.K_GENE], e_w[PredictionModel.WDBQF.K_GENE])
        numpy.testing.assert_array_equal(snp_model.weights[PredictionModel.WDBQF.K_WEIGHT], e_w[PredictionModel.WDBQF.K_WEIGHT])
        numpy.testing.assert_array_equal(snp_model.weights[PredictionModel.WDBQF.K_NON_EFFECT_ALLELE], e_w[PredictionModel.WDBQF.K_NON_EFFECT_ALLELE])
        numpy.testing.assert_array_equal(snp_model.weights[PredictionModel.WDBQF.K_EFFECT_ALLELE], e_w[PredictionModel.WDBQF.K_EFFECT_ALLELE])

    def test_model_manager(self):
        model_manager = PredictionModel.load_model_manager("tests/_td/dbs_2")
        weights = get_weights_in_models("tests/_td/dbs_2")
        models_ = weights[["model", "gene"]].drop_duplicates()
        for t in models_.itertuples():
            compare_model_manager_models_to_weights(self, model_manager, weights, t.model, t.gene)
        #
        self.assertEqual(model_manager.get_genes(), {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'})
        self.assertEqual(model_manager.get_rsids(), set(weights.rsid))
        self.assertEqual(model_manager.get_model_labels(), {"model_sim_1", "model_sim_2"})
        #
        for rsid in set(weights.rsid):
            self.assertEqual(model_manager.snp_keys[rsid], set(weights[weights.rsid == rsid].gene))
        #
        for gene in set(weights.gene):
            w = weights[weights.gene == gene].set_index(["model", "rsid"])[["weight", "eff_allele", "ref_allele"]]
            m = model_manager.get_models(gene).rename(columns={"effect_allele":"eff_allele", "non_effect_allele":"ref_allele"})
            for m_ in set(m.index.get_level_values(0)):
                _compare(self, m.loc[m_], w.loc[m_])

        #
        for gene in set(weights.gene):
            self.assertEqual(model_manager.get_model_labels(gene), set(weights[weights.gene == gene].model))

    def test_optimized_model_manager(self):
        model_manager = PredictionModel.load_model_manager("tests/_td/dbs_2", Klass=PredictionModel._ModelManager)
        weights = get_weights_in_models("tests/_td/dbs_2")

        #
        self.assertEqual(model_manager.get_genes(), {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'})
        self.assertEqual(model_manager.get_rsids(), set(weights.rsid))
        self.assertEqual(model_manager.get_model_labels(), {"model_sim_1", "model_sim_2"})
        #
        for gene in set(weights.gene):
            w = weights[weights.gene == gene]
            m = model_manager.get_models(gene)
            _compare_o(self, m, w)

        #
        for gene in set(weights.gene):
            self.assertEqual(model_manager.get_model_labels(gene), set(weights[weights.gene == gene].model))

    def test_optimized_model_manager_gtex(self):
        model_manager = PredictionModel.load_model_manager("tests/_td/dbs_3", Klass=PredictionModel._ModelManager)
        weights = get_weights_in_models("tests/_td/dbs_3")
        weights.model = weights.model.str.extract("TW_(.*)_0.5", expand=False)

        #
        _assert_optimized_manager(self, model_manager, weights, {"ENSG00000107937.14", "ENSG00000107959.11", "ENSG00000234745.5"})

    def test_optimized_model_manager_gtex_trimmed(self):
        model_manager = PredictionModel.load_model_manager("tests/_td/dbs_3", trim_ensemble_version=True, Klass=PredictionModel._ModelManager)
        weights = get_weights_in_models("tests/_td/dbs_3")
        weights.model = weights.model.str.extract("TW_(.*)_0.5", expand=False)
        weights.gene = weights.gene.str.split(".").str.get(0)

        #
        _assert_optimized_manager(self, model_manager, weights, {"ENSG00000107937", "ENSG00000107959", "ENSG00000234745"})


if __name__ == '__main__':
    unittest.main()
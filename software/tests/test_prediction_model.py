import numpy
import numpy.testing

import unittest
from  metax import PredictionModel

import SampleData

class TestPredictionModelDB(unittest.TestCase):

    def test_load(self):
        t = PredictionModel.ModelDB("tests/_td/dbs/test_1.db")
        extra = t.load_extra()
        self.assertEqual(len(extra), 6)

        e_e = zip(*(SampleData.sample_extra_2()))
        self.assertEqual(extra[PredictionModel.WDBEQF.GENE], e_e[PredictionModel.WDBEQF.GENE])
        self.assertEqual(extra[PredictionModel.WDBEQF.GENE_NAME], e_e[PredictionModel.WDBEQF.GENE_NAME])
        self.assertEqual(extra[PredictionModel.WDBEQF.N_SNP_IN_MODEL], e_e[PredictionModel.WDBEQF.N_SNP_IN_MODEL])
        self.assertEqual(extra[PredictionModel.WDBEQF.PRED_PERF_R2], e_e[PredictionModel.WDBEQF.PRED_PERF_R2])
        self.assertEqual(extra[PredictionModel.WDBEQF.PRED_PERF_PVAL], e_e[PredictionModel.WDBEQF.PRED_PERF_PVAL])
        self.assertEqual(extra[PredictionModel.WDBEQF.PRED_PERF_QVAL], e_e[PredictionModel.WDBEQF.PRED_PERF_QVAL])

        weights = t.load_weights()
        self.assertEqual(len(weights), 5)

        e_w = zip(*(SampleData.sample_weights_2()))
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

if __name__ == '__main__':
    unittest.main()
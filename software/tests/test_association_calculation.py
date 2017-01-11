import numpy
import numpy.testing

import unittest
import SampleData
from metax import PredictionModel
from metax import MatrixManager

from metax.metaxcan import Utilities
from metax.metaxcan import AssociationCalculation

def _gwas():
    gwas = SampleData.dataframe_from_gwas(SampleData.sample_gwas_data_4())
    gwas = gwas.drop(12) #simulate discarding of mismatching alleles
    return gwas

def _prediction_model():
    e = SampleData.dataframe_from_extra(SampleData.sample_extra_1())
    w = SampleData.dataframe_from_weights(SampleData.sample_weights_1())
    p = PredictionModel.Model(w,e)
    return p

def _context():
    gwas = _gwas()
    model = _prediction_model()
    s = SampleData.dataframe_from_covariance(SampleData.sample_covariance_s_1())
    covariance = MatrixManager.MatrixManager(s)
    c = Utilities._build_context(model, covariance, gwas)
    return c

class TestAssociationCalculation(unittest.TestCase):

    def test_intersection_d(self):
        gwas = SampleData.dataframe_from_gwas(SampleData.sample_gwas_data_4())
        p = _prediction_model()
        g, s = Utilities._data_intersection(p, gwas)
        self.assertEqual(set(s), set(['rs1666', 'rs1', 'rs2', 'rs3', 'rs6', 'rs7', 'rs7666', 'rs8', 'rs9','rs100', 'rs101', 'rs102', 'rs202']))
        self.assertEqual(set(g), set(["A","B","C","D"]))

        t_gwas = gwas[3:10]
        g, s = Utilities._data_intersection(p, t_gwas)
        self.assertEqual(set(s), set(['rs3', 'rs6', 'rs7', 'rs7666', 'rs8', 'rs9']))
        self.assertEqual(set(g), set(["B"]))

    def test_build_context(self):
        c = _context()
        r, snps = AssociationCalculation.association("A", c, return_snps=True)
        self.assertEqual(r, ('A', 0.42313735862217716, 0.42845528455235105, 0.10250000000002803, 4, 4, 3))

        r, snps = AssociationCalculation.association("B", c, return_snps=True)
        self.assertEqual(r, ('B', 1.904102672555114, 1.4285714285708686, 0.16333333333323405, 6, 6, 6))

        r, snps = AssociationCalculation.association("C", c, return_snps=True)
        self.assertEqual(r, ('C', 0.089999999999999983, 0.049999999999999989, 0.013333333333320003, 3, 2, 1))

        r, snps = AssociationCalculation.association("D", c, return_snps=True)
        self.assertEqual(r, ('D', numpy.nan, numpy.nan, numpy.nan, 2, numpy.nan, 0))

        r, snps = AssociationCalculation.association("E", c, return_snps=True)
        self.assertEqual(r, ('E', numpy.nan, numpy.nan, numpy.nan, 1, numpy.nan, 0))

    def test_dataframe_from_results(self):
        results = [
            ('A', 0.42313735862217716, 0.42845528455235105, 0.10250000000002803, 4, 4, 3),
            ('B', 1.904102672555114, 1.4285714285708686, 0.16333333333323405, 6, 6, 6),
            ('C', 0.089999999999999983, 0.049999999999999989, 0.013333333333320003, 3, 2, 1)]
        results = zip(*results)
        d = AssociationCalculation.dataframe_from_results(results)
        A = AssociationCalculation.ARF

        numpy.testing.assert_array_equal(d[A.K_GENE], results[A.GENE])
        numpy.testing.assert_array_equal(d[A.K_ZSCORE], results[A.ZSCORE])
        numpy.testing.assert_array_equal(d[A.K_EFFECT_SIZE], results[A.EFFECT_SIZE])
        numpy.testing.assert_array_equal(d[A.K_N_SNPS_IN_MODEL], results[A.N_SNPS_IN_MODEL])
        numpy.testing.assert_array_equal(d[A.K_N_SNPS_IN_COV], results[A.N_SNPS_IN_COV])
        numpy.testing.assert_array_equal(d[A.K_N_SNPS_USED], results[A.N_SNPS_USED])


if __name__ == '__main__':
    unittest.main()
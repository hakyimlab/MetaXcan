import numpy
import numpy.testing

import unittest

from metax import Exceptions
from metax import MatrixManager
D = MatrixManager.GENE_SNP_COVARIANCE_DEFINITION

from . import cov_data
from . import SampleData

class TestMatrixManager(unittest.TestCase):

    def test_invalid_data(self):
        with self.assertRaises(Exceptions.InvalidInputFormat) as ctx:
            MatrixManager.load_matrix_manager("tests/_td/cov/cov.duplicate.txt.gz")

        self.assertTrue("duplicate" in ctx.exception.msg.lower())

        with self.assertRaises(Exceptions.InvalidInputFormat) as ctx:
            MatrixManager.load_matrix_manager("tests/_td/cov/cov.uncontiguous.txt.gz")

        self.assertTrue("contiguous" in ctx.exception.msg.lower())

    def test_from_load(self):
        m = MatrixManager.load_matrix_manager("tests/_td/cov/cov.txt.gz")
        snps, cov = m.get("ENSG00000239789.1")
        self.assertEqual(snps, cov_data.SNPS_ENSG00000239789_1)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_ENSG00000239789_1)

        n = m.n_ids("ENSG00000239789.1")
        self.assertEqual(n, len(cov_data.SNPS_ENSG00000239789_1))

        with self.assertRaises(Exceptions.InvalidArguments) as ctx:
            snps, cov = m.get("ENSG00000183742.8", ["rs7806506", "rs12718973"])

        self.assertTrue("whitelist" in ctx.exception.msg) #?

        whitelist = ["rs3094989", "rs7806506", "rs12536095", "rs10226814"]
        snps, cov = m.get("ENSG00000183742.8", whitelist)
        self.assertEqual(snps, cov_data.SNPS_ENSG00000183742_8_w)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_ENSG00000183742_8_w)

        snps, cov = m.get("ENSG00000004766.11")
        self.assertEqual(snps, cov_data.SNPS_ENSG00000004766_11)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_ENSG00000004766_11)

        n = m.n_ids("ENSG00000004766.11")
        self.assertEqual(n, len(cov_data.COV_ENSG00000004766_11))

    def test_from_data(self):
        s = SampleData.dataframe_from_covariance(SampleData.sample_covariance_s_1())
        m = MatrixManager.MatrixManager(s, D)
        snps, cov = m.get("A")
        self.assertEqual(snps, cov_data.SNPS_A)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_A)

        snps, cov = m.get("B")
        self.assertEqual(snps, cov_data.SNPS_B)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_B)

        snps, cov = m.get("C")
        self.assertEqual(snps, cov_data.SNPS_C)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_C)

        snps, cov = m.get("C", ['rs100', 'rs101', 'rs102'])
        self.assertEqual(snps, cov_data.SNPS_C)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_C)

        with self.assertRaises(Exceptions.InvalidArguments) as ctx:
            snps, cov = m.get("C", ["rs100", "rs12718973"])

        self.assertTrue("whitelist" in ctx.exception.msg)

    def test_flatten(self):
        labels = cov_data.SNPS_ENSG00000183742_8_w
        matrix = cov_data.COV_ENSG00000183742_8_w
        name= "test"

        flat = MatrixManager._flatten_matrix_data([(name, labels, matrix)])
        expected = \
            [('test', 'rs7806506', 'rs7806506', 0.28428631),
             ('test', 'rs7806506', 'rs12536095', -0.01636001),
             ('test', 'rs7806506', 'rs10226814', -0.00157224),
             ('test', 'rs12536095', 'rs12536095', 0.35760734),
             ('test', 'rs12536095', 'rs10226814', 0.00815426),
             ('test', 'rs10226814', 'rs10226814', 0.44923289)]
        numpy.testing.assert_array_equal(flat, expected)

        X = [0,1,3]
        cov = numpy.cov([X])
        flat = MatrixManager._flatten_matrix_data([("a", "b", cov)])

        expected = [('a', 'b', 'b', 2.33333333333333)]

        numpy.testing.assert_array_equal(flat[0][0:3], expected[0][0:3])
        numpy.testing.assert_array_almost_equal(flat[0][3:], expected[0][3:])

if __name__ == '__main__':
    unittest.main()
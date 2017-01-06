import numpy
import numpy.testing

import unittest

from metax import Exceptions
from metax import MatrixManager

import cov_data
import SampleData

class TestMatrixManager(unittest.TestCase):

    def test_invalid(self):
        with self.assertRaises(Exceptions.InvalidInputFormat) as ctx:
            MatrixManager.MatrixManager("tests/_td/cov/cov.duplicate.txt.gz")

        self.assertTrue("duplicate" in ctx.exception.msg)

        with self.assertRaises(Exceptions.InvalidInputFormat) as ctx:
            MatrixManager.MatrixManager("tests/_td/cov/cov.uncontiguous.txt.gz")

        self.assertTrue("contiguous" in ctx.exception.msg)

    def test_load(self):
        m = MatrixManager.MatrixManager("tests/_td/cov/cov.txt.gz")
        snps, cov = m.get("ENSG00000239789.1")
        self.assertEqual(snps, cov_data.SNPS_ENSG00000239789_1)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_ENSG00000239789_1)

        with self.assertRaises(Exceptions.InvalidArguments) as ctx:
            snps, cov = m.get("ENSG00000183742.8", ["rs7806506", "rs12718973"])

        self.assertTrue("whitelist" in ctx.exception.message) #?

        whitelist = ["rs3094989", "rs7806506", "rs12536095", "rs10226814"]
        snps, cov = m.get("ENSG00000183742.8", whitelist)
        self.assertEqual(snps, cov_data.SNPS_ENSG00000183742_8)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_ENSG00000183742_8)

        snps, cov = m.get("ENSG00000004766.11")
        self.assertEqual(snps, cov_data.SNPS_ENSG00000004766_11)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_ENSG00000004766_11)

    def test_from_source(self):
        s = SampleData.dataframe_from_covariance(SampleData.sample_covariance_s_1())
        m = MatrixManager.MatrixManager(s)
        snps, cov = m.get("A")
        self.assertEqual(snps, cov_data.SNPS_A)
        numpy.testing.assert_array_almost_equal(cov, cov_data.COV_A)

if __name__ == '__main__':
    unittest.main()
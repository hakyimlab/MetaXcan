import numpy
import numpy.testing
import pandas

import unittest

from metax import Exceptions
from metax import MatrixManager

import cov_data

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
        self.assertEqual(snps, ['rs7806506', 'rs12536095', 'rs10226814'])
        e_cov = \
        [[0.28428631, -0.01636001, -0.00157224],
         [-0.01636001, 0.35760734, 0.00815426],
         [-0.00157224, 0.00815426, 0.44923289]]
        numpy.testing.assert_array_almost_equal(cov, e_cov)

        snps, cov = m.get("ENSG00000004766.11")
        self.assertEqual(snps, ['rs2285504', 'rs10249649', 'rs10237805', 'rs10235606', 'rs13245529', 'rs7783224', 'rs6979254'])
        e_cov = numpy.matrix(
        [[0.36586853, 0.10533215, -0.00108908, 0.1169794, 0.10533215, 0.09381559, -0.12397725],
        [0.10533215, 0.28500709, 0.00061781, 0.18247091, 0.28500709, 0.0991026, -0.18725892],
        [-0.00108908, 0.00061781, 0.15529136, -0.00195639, -0.00137422, -0.00054256, 0.00362764],
        [0.1169794, 0.18247091, -0.00195639, 0.46566814, 0.18047888, 0.06532518, -0.46544637],
        [0.10533215, 0.28500709, -0.00137422, 0.18047888, 0.28899115, 0.10109463, -0.18725892],
        [0.09381559, 0.0991026, -0.00054256, 0.06532518, 0.10109463, 0.12800488, -0.06759443],
        [-0.12397725, -0.18725892, 0.00362764, -0.46544637, -0.18725892,-0.06759443, 0.47518079]])
        numpy.testing.assert_array_almost_equal(cov, e_cov)

if __name__ == '__main__':
    unittest.main()
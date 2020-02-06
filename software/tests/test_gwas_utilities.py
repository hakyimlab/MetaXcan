import numpy
import numpy.testing
import pandas

import unittest
from metax.gwas import Utilities as GWASUtilities
#
from metax.Constants import SNP
from metax.Constants import EFFECT_ALLELE
from metax.Constants import NON_EFFECT_ALLELE
from metax.Constants import ZSCORE
from metax.Constants import CHROMOSOME
from metax.Constants import POSITION

from . import SampleData

def assert_gwas_1(unit_test, gwas):
    expected_snp = pandas.Series(["rs1666", "rs1", "rs2", "rs3", "rs4", "rs6", "rs7", "rs7666", "rs8", "rs9"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[SNP], expected_snp)

    expected_effect = pandas.Series(["A", "C", "C", "G", "A", "G", "T", "A", "A", "A"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[EFFECT_ALLELE], expected_effect)

    expected_non_effect = pandas.Series(["G", "T", "T", "A", "G", "A", "C", "G", "G", "G"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[NON_EFFECT_ALLELE], expected_non_effect)

    expected_zscore = pandas.Series([0.3, -0.2, 0.5, 1.3, -0.3, 2.9, 4.35, 1.3, 0.09, 0.09], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[ZSCORE], expected_zscore, rtol=0.001)

    expected_chromosome = pandas.Series(["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[CHROMOSOME], expected_chromosome)

    expected_position = pandas.Series([0, 1, 5, 20, 30, 42, 43, 45, 50, 70])
    numpy.testing.assert_array_equal(gwas[POSITION], expected_position)

def assert_gwas_2(unit_test, gwas):
    expected_snp = pandas.Series(["rsC", "rs1666", "rs1", "rs2",  "rs4", "rsB", "rsA", "rs7666", "rs8", "rs9"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[SNP], expected_snp)

    expected_effect = pandas.Series(["T", "A", "C", "C", "A", "G", "G", "A", "A", "A"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[EFFECT_ALLELE], expected_effect)

    expected_non_effect = pandas.Series(["C", "G", "T", "T", "G", "A", "A", "G", "G", "G"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[NON_EFFECT_ALLELE], expected_non_effect)

    expected_zscore = pandas.Series([4.35, 0.3, -0.2, 1.3, -0.3, 2.9, 1.3, 1.3, 0.09, 0.09], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[ZSCORE], expected_zscore, rtol=0.001)

    expected_chromosome = pandas.Series(["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[CHROMOSOME], expected_chromosome)

    expected_position = pandas.Series([None, None, None, None, None, None, None, None, None, None])
    numpy.testing.assert_array_equal(gwas[POSITION], expected_position)

def assert_gwas_1_e(unit_test, gwas):
    assert_gwas_1(unit_test, gwas)

    numpy.testing.assert_array_equal(gwas["number"], numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))

    numpy.testing.assert_array_equal(gwas["character"], numpy.array(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]))

class TestGWASUtilities(unittest.TestCase):

    def test_gwas_from_data(self):
        gwas = GWASUtilities.gwas_from_data(SampleData.sample_gwas_data_1())
        assert_gwas_1(self, gwas)

        gwas = GWASUtilities.gwas_from_data(SampleData.sample_gwas_data_2())
        assert_gwas_2(self, gwas)

        gwas = GWASUtilities.gwas_from_data(SampleData.sample_gwas_data_1_e(), extra_columns=[("number",6), ("character",7)])
        assert_gwas_1_e(self, gwas)

if __name__ == '__main__':
    unittest.main()
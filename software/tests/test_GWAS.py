import sys

import numpy
import numpy.testing
import pandas

import unittest

from metax.gwas import GWAS
from metax.gwas import Utilities as GWASUtilities

from metax.Constants import SNP
from metax.Constants import EFFECT_ALLELE
from metax.Constants import NON_EFFECT_ALLELE
from metax.Constants import ZSCORE
from metax.Constants import CHROMOSOME
from metax.Constants import POSITION

import SampleData

def assert_basic_gwas(unit_test, gwas):
    expected_snp = pandas.Series(["rs940550", "rs6650104", "rs6594028", "rs9701055", "rs7417504", "rs12082473", "rs3094315", "rs3131971", "rs61770173", ], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[SNP], expected_snp)

    expected_effect = pandas.Series(["C", "T", "T", "A", "T", "A", "A", "T", "A"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[EFFECT_ALLELE], expected_effect)

    expected_non_effect = pandas.Series(["G", "C", "C", "T", "C", "G", "G", "C", "C"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[NON_EFFECT_ALLELE], expected_non_effect)

    expected_zscore = pandas.Series([-1.254557, 0.974874, 1.024923, -0.652800, -0.194823, -0.268321, 0.226338, 0.666656, -0.232505], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[ZSCORE], expected_zscore, rtol=0.001)

def assert_gwas_with_chromosome(unit_test, gwas):
    assert_basic_gwas(unit_test, gwas)

    expected_chromosome = pandas.Series(["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr22", "chr22"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[CHROMOSOME], expected_chromosome)


def assert_full_gwas(unit_test, gwas):
    assert_gwas_with_chromosome(unit_test, gwas)

    expected_position = pandas.Series([729679, 731718, 734349, 736289, 751756, 752566, 752721, 752894, 753405])
    numpy.testing.assert_array_equal(gwas[POSITION], expected_position)


def assert_gwas_extracted_from_data_3(unit_test, gwas):
    expected_snp = pandas.Series(["rs3", "rs6", "rs7"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[SNP], expected_snp)

    expected_effect = pandas.Series(["G", "G", "T"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[EFFECT_ALLELE], expected_effect)

    expected_non_effect = pandas.Series(["A", "A", "C"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[NON_EFFECT_ALLELE], expected_non_effect)

    expected_zscore = pandas.Series([1.3, 2.9, 4.35], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[ZSCORE], expected_zscore, rtol=0.001)

    expected_chromosome = pandas.Series(["chr1", "chr1", "chr1"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[CHROMOSOME], expected_chromosome)

class TestGWAS(unittest.TestCase):

    def test_load_gwas(self):
        #full format
        gwas_format = {
            "column_snp":"snpid",
            "column_non_effect_allele":"a2",
            "column_effect_allele":"a1",
            "column_or":"or",
            "column_se":"se",
            "column_chromosome":"hg19chrc",
            "column_position":"bp"
        }

        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_full_gwas(self, gwas)

        # missing fields
        gwas_format = {
            "column_snp":"snpid",
            "column_non_effect_allele":"a2",
            "column_effect_allele":"a1",
            "column_or":"or",
            "column_se":"se",
        }

        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_basic_gwas(self, gwas)

    def test_extract(self):
        gwas = GWASUtilities.gwas_from_data(SampleData.sample_gwas_data_3())
        g = GWAS.extract(gwas, ["rs3", "rs6", "rs7"])
        assert_gwas_extracted_from_data_3(self, g)

if __name__ == "__main__":
    unittest.main()
import numpy
import numpy.testing
import pandas

import unittest

from metax.gwas import GWAS
from metax.gwas import Utilities as GWASUtilities

from metax.Constants import SNP
from metax.Constants import EFFECT_ALLELE
from metax.Constants import NON_EFFECT_ALLELE
from metax.Constants import PVALUE
from metax.Constants import BETA
from metax.Constants import SE
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


def assert_gwas_with_chromosome(unit_test, gwas):
    assert_basic_gwas(unit_test, gwas)

    expected_chromosome = pandas.Series(["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr22", "chr22"], dtype=numpy.str)
    numpy.testing.assert_array_equal(gwas[CHROMOSOME], expected_chromosome)

def assert_extended_gwas(unit_test, gwas):
    assert_gwas_with_chromosome(unit_test, gwas)

    expected_position = pandas.Series([729679, 731718, 734349, 736289, 751756, 752566, 752721, 752894, 753405])
    numpy.testing.assert_array_equal(gwas[POSITION], expected_position)

def assert_gwas_zscore_fbse(unit_test, gwas):
    assert_extended_gwas(unit_test, gwas)

    expected_zscore = pandas.Series( [-1.254557, 0.974874, 1.024923, -0.652800, -0.194823, -0.268321, 0.226338, 0.666656, -0.232505], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[ZSCORE], expected_zscore, rtol=0.001)

    expected_beta = pandas.Series( [-0.0217038334437866, 0.0193025022544974, 0.0204984635773248, -0.0125990355765152, -0.00319509889654086, -0.00399798128729837, 0.00330453400830047, 0.0099998345783334, -0.00369682484428976], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[BETA], expected_beta, rtol=0.001)

    expected_se = pandas.Series( [0.0173, 0.0198, 0.02, 0.0193, 0.0164, 0.0149, 0.0146, 0.015, 0.0159], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[SE], expected_se, rtol=0.001)

def assert_gwas_zscore_pb(unit_test, gwas):
    assert_extended_gwas(unit_test, gwas)

    expected_zscore = pandas.Series( [-1.258254,  0.974517,  1.02471 , -0.653863, -0.19793 , -0.270208, 0.223817,  0.664297, -0.229989], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[ZSCORE], expected_zscore, rtol=0.001)

    expected_p = pandas.Series( [0.2083, 0.3298, 0.3055, 0.5132, 0.8431, 0.7870, 0.8229, 0.5065, 0.8181], dtype=numpy.float32)
    numpy.testing.assert_allclose(gwas[PVALUE], expected_p, rtol=0.001)

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
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format, strict=False)
        assert_basic_gwas(self, gwas)

        #full format, OR+SE (which is like beta+se)
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_or":"OR",
            "column_se":"SE",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_gwas_zscore_fbse(self, gwas)

        # full format, beta+SE
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_beta":"BETA",
            "column_se":"SE",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_gwas_zscore_fbse(self, gwas)

        # full format, pvalue+beta
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_beta":"BETA",
            "column_pvalue":"P",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_gwas_zscore_pb(self, gwas)

        # full format, pvalue+beta_sign
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_beta_sign":"BETA_SIGN",
            "column_pvalue":"P",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_gwas_zscore_pb(self, gwas)

        # full format, pvalue+or
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_beta_sign":"BETA_SIGN",
            "column_pvalue":"P",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2.gwas.results.txt.gz", gwas_format)
        assert_gwas_zscore_pb(self, gwas)

    def test_gwas_from_source(self):
        #full format, OR+SE (which is like beta+se)
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_or":"OR",
            "column_se":"SE",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }

        source = GWASUtilities.gwas_filtered_source("tests/_td/GWAS/scz2.gwas.results.txt.gz")
        gwas = GWAS.load_gwas(source, gwas_format)
        assert_gwas_zscore_fbse(self, gwas)

        source = GWASUtilities.gwas_filtered_source("tests/_td/GWAS/scz2.gwas.results.txt.gz", snps={"rs940550", "rs6650104", "rs61770173"}, snp_column_name="SNPID")
        gwas = GWAS.load_gwas(source, gwas_format)

        numpy.testing.assert_array_equal(gwas[SNP], pandas.Series(["rs940550", "rs6650104", "rs61770173", ], dtype=numpy.str))
        numpy.testing.assert_array_equal(gwas[EFFECT_ALLELE], pandas.Series(["C", "T",  "A"], dtype=numpy.str))
        numpy.testing.assert_array_equal(gwas[NON_EFFECT_ALLELE], pandas.Series(["G", "C", "C"], dtype=numpy.str))
        numpy.testing.assert_array_equal(gwas[CHROMOSOME], pandas.Series(["chr1", "chr1",  "chr22"], dtype=numpy.str))
        numpy.testing.assert_allclose(gwas[ZSCORE], pandas.Series([-1.254557, 0.974874, -0.232505],dtype=numpy.float32), rtol=0.001)
        numpy.testing.assert_allclose(gwas[BETA], pandas.Series([-0.0217038334437866, 0.0193025022544974, -0.00369682484428976], dtype=numpy.float32), rtol=0.001)
        numpy.testing.assert_allclose(gwas[SE], pandas.Series([0.0173, 0.0198,  0.0159], dtype=numpy.float32), rtol=0.001)


    def test_extract(self):
        gwas = GWASUtilities.gwas_from_data(SampleData.sample_gwas_data_3())
        g = GWAS.extract(gwas, ["rs3", "rs6", "rs7"])
        assert_gwas_extracted_from_data_3(self, g)


if __name__ == "__main__":
    unittest.main()
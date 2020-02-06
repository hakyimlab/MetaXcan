import numpy
import numpy.testing
import pandas

import unittest

from metax import Exceptions

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

from . import SampleData
from . import scz2_sample

def assert_basic_gwas(unit_test, gwas):
    numpy.testing.assert_array_equal(gwas[SNP], scz2_sample.expected_snp)
    numpy.testing.assert_array_equal(gwas[EFFECT_ALLELE], scz2_sample.expected_effect)
    numpy.testing.assert_array_equal(gwas[NON_EFFECT_ALLELE], scz2_sample.expected_non_effect)


def assert_gwas_with_chromosome(unit_test, gwas):
    assert_basic_gwas(unit_test, gwas)
    numpy.testing.assert_array_equal(gwas[CHROMOSOME], scz2_sample.expected_chromosome)

def assert_extended_gwas(unit_test, gwas):
    assert_gwas_with_chromosome(unit_test, gwas)
    numpy.testing.assert_array_equal(gwas[POSITION], scz2_sample.expected_position)

def assert_gwas_zscore_fbse(unit_test, gwas):
    assert_extended_gwas(unit_test, gwas)
    numpy.testing.assert_allclose(gwas[ZSCORE], scz2_sample.expected_zscore_1, rtol=0.001)
    numpy.testing.assert_allclose(gwas[BETA], scz2_sample.expected_beta, rtol=0.001)
    numpy.testing.assert_allclose(gwas[SE], scz2_sample.expected_se, rtol=0.001)

def assert_gwas_zscore_pb(unit_test, gwas):
    assert_extended_gwas(unit_test, gwas)
    numpy.testing.assert_allclose(gwas[ZSCORE], scz2_sample.expected_zscore_2, rtol=0.001)
    numpy.testing.assert_allclose(gwas[PVALUE], scz2_sample.expected_p, rtol=0.001)

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

def _add_basic_to_format(format):
    format[GWAS.COLUMN_SNP] = "_snp"
    format[GWAS.COLUMN_EFFECT_ALLELE] = "_effect_allele"
    format[GWAS.COLUMN_NON_EFFECT_ALLELE] = "_non_effect_allele"

def _add_extra_to_format(format):
    format[GWAS.COLUMN_ZSCORE] = "_zscore"
    format[GWAS.COLUMN_BETA] = "_beta"
    format[GWAS.COLUMN_BETA_SIGN] = "_beta_sign"
    format[GWAS.COLUMN_OR] = "_or"
    format[GWAS.COLUMN_SE] = "_se"
    format[GWAS.COLUMN_PVALUE] = "_pvalue"

class TestGWAS(unittest.TestCase):

    def test_load_gwas(self):
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format, strict=False)
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
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format)
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
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format)
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
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format)
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
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format)
        assert_gwas_zscore_pb(self, gwas)

        # full format, pvalue+or
        gwas_format = {
            "column_snp":"SNPID",
            "column_non_effect_allele":"A2",
            "column_effect_allele":"A1",
            "column_or":"OR",
            "column_pvalue":"P",
            "column_chromosome":"HG19CHRC",
            "column_position":"BP"
        }
        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format)
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

        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format, force_special_handling=True)
        assert_gwas_zscore_fbse(self, gwas)

        gwas = GWAS.load_gwas("tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz", gwas_format, snps={"rs940550", "rs6650104", "rs61770173"})

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

    def test_format(self):
        format = {}
        self.assertIsNone(GWAS._f_snp(format))
        self.assertIsNone(GWAS._f_effect_allele_column(format))
        self.assertIsNone(GWAS._f_non_effect_allele_column(format))
        self.assertIsNone(GWAS._f_pvalue(format))
        self.assertIsNone(GWAS._f_zscore(format))
        self.assertIsNone(GWAS._f_beta(format))
        self.assertIsNone(GWAS._f_beta_sign(format))
        self.assertIsNone(GWAS._f_or(format))
        self.assertIsNone(GWAS._f_se(format))

        _add_basic_to_format(format)
        _add_extra_to_format(format)

        self.assertEqual(GWAS._f_snp(format), "_snp")
        self.assertEqual(GWAS._f_effect_allele_column(format), "_effect_allele")
        self.assertEqual(GWAS._f_non_effect_allele_column(format), "_non_effect_allele")
        self.assertEqual(GWAS._f_pvalue(format), "_pvalue")
        self.assertEqual(GWAS._f_zscore(format), "_zscore")
        self.assertEqual(GWAS._f_beta(format), "_beta")
        self.assertEqual(GWAS._f_beta_sign(format), "_beta_sign")
        self.assertEqual(GWAS._f_or(format), "_or")
        self.assertEqual(GWAS._f_se(format), "_se")

    def test_format_validation(self):
        format = {}
        with self.assertRaises(Exceptions.InvalidArguments): GWAS.validate_format_basic(format)
        with self.assertRaises(Exceptions.InvalidArguments): GWAS.validate_format_for_strict(format)

        _add_basic_to_format(format)
        GWAS.validate_format_basic(format)
        with self.assertRaises(Exceptions.InvalidArguments): GWAS.validate_format_for_strict(format)

        format[GWAS.COLUMN_PVALUE] = "_p"
        format[GWAS.COLUMN_SE] = "_se"
        with self.assertRaises(Exceptions.InvalidArguments): GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_ZSCORE] = "_z"
        GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_PVALUE] = "_p"
        format[GWAS.COLUMN_BETA] = "_beta"
        GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_PVALUE] = "_p"
        format[GWAS.COLUMN_BETA] = "_beta_sign"
        GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_PVALUE] = "_p"
        format[GWAS.COLUMN_OR] = "_or"
        GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_SE] = "se"
        with self.assertRaises(Exceptions.InvalidArguments): GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_SE] = "se"
        format[GWAS.COLUMN_OR] = "or"
        GWAS.validate_format_for_strict(format)

        #
        format = {}
        _add_basic_to_format(format)
        format[GWAS.COLUMN_SE] = "se"
        format[GWAS.COLUMN_OR] = "beta"
        GWAS.validate_format_for_strict(format)

if __name__ == "__main__":
    unittest.main()
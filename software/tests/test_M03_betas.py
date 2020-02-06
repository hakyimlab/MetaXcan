#!/usr/bin/env python
import unittest
import shutil
import os

import logging
import numpy
import numpy.testing
import pandas

from unittest.mock import patch

from metax.Constants import SNP
from metax.Constants import BETA
from metax.Constants import ZSCORE

from metax import Exceptions

from M03_betas import run
from . import scz2_sample

class DummyArgs(object):
    def __init__(self):
        self.snp_column = None
        self.effect_allele_column = None
        self.non_effect_allele_column = None
        self.chromosome_column = None
        self.position_column = None
        self.beta_column = None
        self.beta_sign_column = None
        self.se_column = None
        self.or_column = None
        self.zscore_column = None
        self.freq_column = None
        self.pvalue_column = None
        self.gwas_folder = None
        self.gwas_file_pattern = None
        self.output_folder = None
        self.separator = None
        self.skip_until_header = None
        self.verbosity = logging.CRITICAL
        self.throw = True
        self.model_db_path = None
        self.handle_empty_columns = False
        self.input_pvalue_fix = 1e-30
        self.gwas_file = None
        self.model_db_snp_key = None
        self.keep_non_rsid = None
        self.output = None
        self.snp_map_file = None

def base_args(folder="tests/_td/GWAS/scz2", file=None):
    args = DummyArgs()
    args.gwas_folder = folder
    args.gwas_file = file
    args.snp_column = "SNPID"
    args.effect_allele_column = "A1"
    args.non_effect_allele_column = "A2"
    return args

def assert_beta_p(unit_test, r, filter=None):
    numpy.testing.assert_array_equal(r[SNP], scz2_sample.expected_snp)
    numpy.testing.assert_allclose(r[ZSCORE], scz2_sample.expected_zscore_2, rtol=0.001)

def assert_beta_pb(unit_test, r):
    assert_beta_p(unit_test, r)
    numpy.testing.assert_allclose(r[BETA], scz2_sample.expected_beta, rtol=0.001)

def assert_beta_bse(unit_test, r):
    numpy.testing.assert_array_equal(r[SNP], scz2_sample.expected_snp)
    numpy.testing.assert_allclose(r[ZSCORE], scz2_sample.expected_zscore_1, rtol=0.001)
    numpy.testing.assert_allclose(r[BETA], scz2_sample.expected_beta, rtol=0.001)

def assert_beta_pb_fix(unit_test, r):
    numpy.testing.assert_array_equal(r[SNP], scz2_sample.expected_snp)
    numpy.testing.assert_allclose(r[ZSCORE], scz2_sample.expected_zscore_fix, rtol=0.001)
    numpy.testing.assert_allclose(r[BETA], scz2_sample.expected_beta_fix, rtol=0.001)

def assert_model_beta_pb(unit_test, r):
    expected_snp = pandas.concat([scz2_sample.expected_snp[0:3],scz2_sample.expected_snp[-2:]]).values
    numpy.testing.assert_array_equal(r[SNP], expected_snp)

    expected_zscore = pandas.concat([scz2_sample.expected_zscore_2[0:3],scz2_sample.expected_zscore_2[-2:]]).values
    expected_zscore[2] = -expected_zscore[2]
    numpy.testing.assert_allclose(r[ZSCORE], expected_zscore, rtol=0.001)

    expected_beta = pandas.concat([scz2_sample.expected_beta[0:3],scz2_sample.expected_beta[-2:]]).values
    expected_beta[2] = -expected_beta[2]
    numpy.testing.assert_allclose(r[BETA], expected_beta, rtol=0.001)

class TestM03(unittest.TestCase):

    @patch('metax.gwas.GWAS.validate_format_for_strict')
    @patch('metax.gwas.GWAS.validate_format_basic')
    def testWrongArguments(self, patch_validate_basic, patch_validate_strict):
        args = DummyArgs()
        with self.assertRaises(Exception) as c: run(args)

        args.gwas_folder = "tests/_td/GWAS/scz2"
        patch_validate_basic.side_effect = RuntimeError("k")
        with self.assertRaises(RuntimeError) as c: run(args)

        patch_validate_basic.side_effect = None
        patch_validate_strict.side_effect = RuntimeError()
        with self.assertRaises(RuntimeError) as c: run(args)

        patch_validate_basic.side_effect = None
        patch_validate_strict.side_effect = None

    def test_fail_incompatible_arguments(self):
        args = base_args("tests/_td/GWAS/scz2", "tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz")
        with self.assertRaises(Exceptions.InvalidArguments) as c: run(args)

    def test_run_folder(self):
        self.run_folder_or_file("tests/_td/GWAS/scz2", None)

    def test_run_file(self):
        self.run_folder_or_file(None, "tests/_td/GWAS/scz2/scz2.gwas.results.txt.gz")

    def run_folder_or_file(self, folder, file):
        args = base_args(folder, file)
        args.pvalue_column = "P"
        args.or_column = "OR"
        r = run(args)
        assert_beta_pb(self, r)

        args = base_args(folder, file)
        args.pvalue_column = "P"
        args.beta_column = "BETA"
        r = run(args)
        assert_beta_pb(self, r)

        args = base_args(folder, file)
        args.pvalue_column = "P"
        args.beta_sign_column = "BETA_SIGN"
        r = run(args)
        assert_beta_p(self, r)

        args = base_args(folder, file)
        args.beta_column = "BETA"
        args.se_column = "SE"
        r = run(args)
        assert_beta_bse(self, r)

        args = base_args(folder, file)
        args.beta_column = "BETA"
        args.se_column = "SE"
        r = run(args)
        assert_beta_bse(self, r)

        args = base_args(folder, file)
        args.or_column = "OR"
        args.se_column = "SE"
        r = run(args)
        assert_beta_bse(self, r)

        #Should fail
        args = base_args(folder, file)
        args.or_column = "BETA"
        args.se_column = "SE"
        with self.assertRaises(Exception): M03_betas.run(args)

    def test_run_skip(self):
        args = base_args("tests/_td/GWAS/scz2c")
        args.pvalue_column = "P"
        args.or_column = "OR"
        args.skip_until_header = "\t".join(["HG19CHRC", "SNPID", "A1", "A2", "BP", "INFO", "OR", "SE", "P", "NGT", "BETA", "ZSCORE", "BETA_SIGN"])
        r = run(args)
        assert_beta_pb(self, r)

    def test_run_split(self):
        args = base_args("tests/_td/GWAS/scz2b")
        args.gwas_file_pattern = ".*gz"
        args.pvalue_column = "P"
        args.or_column = "OR"
        r = run(args)
        assert_beta_pb(self, r)

    def test_with_model(self):
        args = base_args()
        args.pvalue_column = "P"
        args.or_column = "OR"
        args.model_db_path = "tests/_td/dbs/test_3.db"
        r = run(args)

        assert_model_beta_pb(self, r)

    def test_split_with_model(self):
        op = ".kk_test"
        if os.path.exists(op): shutil.rmtree(op)
        args = base_args("tests/_td/GWAS/scz2b")
        args.gwas_file_pattern = ".*gz"
        args.pvalue_column = "P"
        args.or_column = "OR"
        args.model_db_path = "tests/_td/dbs/test_3.db"
        r = run(args)

        assert_model_beta_pb(self, r)

    def test_run_to_file(self):
        op = ".kk_test"
        if os.path.exists(op): shutil.rmtree(op)
        args = base_args()
        args.pvalue_column = "P"
        args.or_column = "OR"
        args.output_folder = op
        run(args)

        n = os.listdir(op)[0]
        n = os.path.join(op,n)
        r = pandas.read_table(n)

        assert_beta_pb(self, r)
        shutil.rmtree(op)

    def test_run_to_files(self):
        op = ".kk_test"
        if os.path.exists(op): shutil.rmtree(op)
        args = base_args("tests/_td/GWAS/scz2b")
        args.gwas_file_pattern = ".*gz"
        args.pvalue_column = "P"
        args.or_column = "OR"
        args.output_folder = op
        run(args)

        r = pandas.DataFrame()
        for n in sorted(os.listdir(op)):
            n = os.path.join(op,n)
            d = pandas.read_table(n)
            r = pandas.concat([r,d])

        assert_beta_pb(self, r)
        shutil.rmtree(op)

    def test_fix_pvalue(self):
        op = ".kk_test"
        if os.path.exists(op): shutil.rmtree(op)
        args = base_args("tests/_td/GWAS/scz2d")
        args.gwas_file_pattern = ".*gz"
        args.pvalue_column = "P"
        args.or_column = "OR"
        r = run(args)
        assert_beta_pb_fix(self, r)

    def test_align(self):
        pass

if __name__ == "__main__":
    unittest.main()
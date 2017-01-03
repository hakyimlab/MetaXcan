#!/usr/bin/env python
import unittest
import sys
import shutil
import os
import re
import gzip
import logging

from mock import patch
from mock import Mock
from mock import call

from metax import Exceptions

import M03_betas
import test_gwas


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
        self.verbosity = logging.ERROR
        self.throw = True
        self.model_db_path = None

class TestM03(unittest.TestCase):

    @patch('metax.gwas.GWAS.validate_format_for_strict')
    @patch('metax.gwas.GWAS.validate_format_basic')
    def testWrongArguments(self, patch_validate_basic, patch_validate_strict):
        args = DummyArgs()
        with self.assertRaises(Exception) as c: M03_betas.run(args)

        args.gwas_folder = "tests/_td/GWAS/scz2"
        patch_validate_basic.side_effect = Exception()
        with self.assertRaises(Exception) as c: M03_betas.run(args)

        patch_validate_basic.side_effect = None
        patch_validate_strict.side_effect = Exception()
        with self.assertRaises(Exception) as c: M03_betas.run(args)

if __name__ == "__main__":
    unittest.main()
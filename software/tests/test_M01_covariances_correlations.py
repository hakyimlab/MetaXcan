#!/usr/bin/env python
import unittest
import sys
import shutil
import os
import re
import gzip

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.Formats as Formats

from M01_covariances_correlations import ProcessWeightDB

class Dummy(object):
    pass

def buildDummyArgs(root):
    dummy = Dummy()
    dummy.verbosity = 10
    dummy.weight_db = os.path.join(root, "test.db")
    dummy.input_folder = os.path.join(root, "filtered_dosage")
    dummy.covariance_output = os.path.join(root, "covariance/cov.txt.gz")
    dummy.correlation_output = None
    dummy.input_format = Formats.PrediXcan
    dummy.min_maf_filter = None
    dummy.max_maf_filter = None
    dummy.max_snps_in_gene = None
    dummy.delimiter = " "
    return dummy

def setupDataForArgs(args, root):
    if os.path.exists(root):
        shutil.rmtree(root)
    shutil.copytree("tests/_td/filtered_dosage", os.path.join(root,"filtered_dosage"))
    shutil.copy("tests/_td/test.db", root)

def cleanUpDataForArgs(root):
    shutil.rmtree(root)

class TestM01(unittest.TestCase):
    def testProcessPrerequisitesnoArgConstructor(self):
        with self.assertRaises(AttributeError):
            dummy = Dummy()
            p = ProcessWeightDB(dummy)

    def testProcessPrerequisitesConstructor(self):
        dummy = buildDummyArgs("_test")
        p = ProcessWeightDB(dummy)
        self.assertEqual(p.weight_db, "test.db")
        self.assertEqual(p.db_path, "_test/test.db")
        self.assertEqual(p.data_folder, "_test/filtered_dosage")
        self.assertEqual(p.correlation_output, None)
        self.assertEqual(p.covariance_output, "_test/covariance/cov.txt.gz")
        self.assertEqual(p.input_format, Formats.PrediXcan)
        self.assertEqual(p.min_maf_filter, None)
        self.assertEqual(p.max_maf_filter, None)

    def testProcessPrerequisitesConstructorDefaultCovariance(self):
        dummy = buildDummyArgs("_test")
        dummy.covariance_output = None
        p = ProcessWeightDB(dummy)
        self.assertEqual(p.covariance_output, "intermediate/cov/test.cov.txt.gz")

    def testProcessWeightDBRun(self):
        dummy = buildDummyArgs("_test")
        setupDataForArgs(dummy, "_test")
        p = ProcessWeightDB(dummy)

        try:
            p.run()
        except:
            self.assertEqual(False, True, "Prerequisites should have run without error")

        with gzip.open(p.covariance_output) as f:
            expected_lines = ["GENE RSID1 RSID2 VALUE",
                            "A rs1 rs1 1.0",
                            "A rs1 rs2 0.0",
                            "A rs1 rs3 -0.3333333333333333",
                            "A rs2 rs2 0.0",
                            "A rs2 rs3 0.0",
                            "A rs3 rs3 0.3333333333333333",
                            "B rs4 rs4 0.0",
                            "B rs4 rs5 0.0",
                            "B rs5 rs5 0.3333333333333333",
                            "C rs6 rs6 1.0",
                            "D rs1 rs1 1.0"]
            for i,expected_line in enumerate(expected_lines):
                actual_line = f.readline().decode().strip()
                self.assertEqual(actual_line, expected_line)
        cleanUpDataForArgs("_test")

#!/usr/bin/env python
import unittest
import sys
import shutil
import os
import re

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

from metax.DataSet import DataSet
from metax.DataSet import DataSetFileUtilities
from metax.DataSet import DataSetCollection
import metax.Formats as Formats

from M00_prerequisites import ProcessPrerequisites

class Dummy(object):
    pass

def buildDummyArgs(root):
    dummy = Dummy()
    dummy.verbosity = 10
    dummy.dosage_folder = os.path.join(root, "dosage_set_1")
    dummy.snp_list = os.path.join(root, "snp.list.gz")
    dummy.output_folder = os.path.join(root, "intermediate/filtered")
    dummy.file_pattern = "set_(.*)"
    dummy.population_filters = ["HERO"]
    dummy.individual_filters = []
    dummy.input_format = Formats.IMPUTE
    dummy.output_format = Formats.PrediXcan
    return dummy

def setupDataForArgs(args, root):
    if os.path.exists(root):
        shutil.rmtree(root)
    shutil.copytree("tests/_td", root)

def cleanUpDataForArgs(root):
    shutil.rmtree(root)

class TestM00(unittest.TestCase):

    def testProcessPrerequisitesnoArgConstructor(self):
        with self.assertRaises(AttributeError):
            dummy = Dummy()
            p = ProcessPrerequisites(dummy)

    def testProcessPrerequisitesConstructor(self):
        dummy = buildDummyArgs("_test")
        setupDataForArgs(dummy, "_test")
        p = ProcessPrerequisites(dummy)
        self.assertEqual(p.dosage_folder, "_test/dosage_set_1")
        self.assertEqual(p.snp_list, "_test/snp.list.gz")
        self.assertEqual(p.output_folder, "_test/intermediate/filtered")
        self.assertEqual(p.population_filters, ["HERO"])
        self.assertEqual(p.individual_filters, [])
        self.assertEqual(p.chromosome_in_name_regex,re.compile("set_(.*)"))
        self.assertEqual(p.samples_input, "_test/dosage_set_1/set.sample")
        self.assertEquals(p.samples_output, "_test/intermediate/filtered/set.sample")
        cleanUpDataForArgs("_test")

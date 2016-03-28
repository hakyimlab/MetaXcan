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

import metax.ThousandGenomesUtilities as ThousandGenomesUtilities

def setupData(root):
    if os.path.exists(root):
        shutil.rmtree(root)
    shutil.copytree("tests/_td/dosage_set_1", os.path.join(root,"dosage_set_1"))
    shutil.copy("tests/_td/snp.txt.gz", root)

def cleanUpData(root):
    shutil.rmtree(root)

class TestThousandGenomesUtilities(unittest.TestCase):
    def testLegendLoader(self):
        class DummyCalback:
            def __init__(self):
                self.lines = []

            def __call__(self, index, row):
                self.lines.append((index, row))

        callback = DummyCalback()
        loader = ThousandGenomesUtilities.LEGENDLoader("tests/_td/dosage_set_1", "set_chr1")
        loader.iterateOverFileLegends(callback)
        self.assertEqual(6, len(callback.lines))
        expected = [
            (0, ["1:10177:A:AC", "10177", "A", "AC", "Biallelic_INDEL", "0.490922844175492", "0.360230547550432", "0.336309523809524", "0.405566600397614", "0.494887525562372", "0.425319488817891"]),
            (1, ["rs1:1:A:T", "10505", "A", "T", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (2, ["1:12:C:G", "10506", "C", "G", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (3, ["rs2:2:G:A", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (4, ["rs3:3:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (5, ["rs4:4:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"])
        ]

        for i, line in enumerate(callback.lines):
            e = expected[i]
            self.assertEqual(e, line)
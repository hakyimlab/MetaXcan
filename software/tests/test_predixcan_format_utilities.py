import unittest
import sys
import shutil
import os
import gzip
import re

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.WeightDBUtilities as WeightDBUtilities

class TestPrediXcanFormatUtilities(unittest.TestCase):
    def testPrediXcanLoader(self):
        weight_db = WeightDBUtilities.WeightDBEntryLogic("tests/_td/test.db")
        loader = PrediXcanFormatUtilities.PrediXcanFormatDosageLoader("tests/_td/filtered_dosage/chr1.dosage.gz", weight_db)
        snps, snps_by_rsid = loader.load()

        self.assertEqual(sorted(snps_by_rsid.keys()), ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"])

        self.assertEqual(snps[0].ref_allele, "T")
        self.assertEqual(snps[0].eff_allele, "C")
        self.assertEqual(snps[0].position, 1)
        self.assertEqual(snps[0].data, [0, 0, 0, 2])

        self.assertEqual(snps[1].ref_allele, "A")
        self.assertEqual(snps[1].eff_allele, "G")
        self.assertEqual(snps[1].position, 2)
        self.assertEqual(snps[1].data, [1, 1, 1, 1])

        self.assertEqual(snps[2].ref_allele, "G")
        self.assertEqual(snps[2].eff_allele, "A")
        self.assertEqual(snps[2].position, 3)
        self.assertEqual(snps[2].data, [0, 1, 0, 1])

        self.assertEqual(snps[3].ref_allele, "T")
        self.assertEqual(snps[3].eff_allele, "C")
        self.assertEqual(snps[3].position, 4)
        self.assertEqual(snps[3].data, [0, 0, 0, 0])

        self.assertEqual(snps[4].ref_allele, "C")
        self.assertEqual(snps[4].eff_allele, "T")
        self.assertEqual(snps[4].position, 5)
        self.assertEqual(snps[4].data, [0, 0, 1, 1])

        self.assertEqual(snps[5].ref_allele, "C")
        self.assertEqual(snps[5].eff_allele, "T")
        self.assertEqual(snps[5].position, 6)
        self.assertEqual(snps[5].data, [0, 0, 0, 2])

        # 7th row doesn't get loaded because it doesn't belong to weight db

if __name__ == '__main__':
    unittest.main()
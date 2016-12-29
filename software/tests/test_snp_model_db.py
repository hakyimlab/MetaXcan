import sys

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest
from  metax import SNPModelDB

class TestTranscriptomeDB(unittest.TestCase):

    def test_load_dosage(self):
        t = SNPModelDB.SNPModelDB("tests/_td/dbs/test_1.db")
        extra = t.load_extra()
        self.assertEqual(len(extra), 6)

        self.assertEqual(extra[SNPModelDB.WDBEQF.GENE], (u'A', u'B', u'C', u'D'))
        self.assertEqual(extra[SNPModelDB.WDBEQF.GENE_NAME], (u'gene1', u'gene2', u'gene3', u'gene4'))
        self.assertEqual(extra[SNPModelDB.WDBEQF.N_SNP_IN_MODEL], (3, 2, 1, 1))
        self.assertEqual(extra[SNPModelDB.WDBEQF.PRED_PERF_R2], (0.9, 0.8, 0.7, 0.6))
        self.assertEqual(extra[SNPModelDB.WDBEQF.PRED_PERF_PVAL], (0.09, 0.08, 0.07, 0.06))
        self.assertEqual(extra[SNPModelDB.WDBEQF.PRED_PERF_QVAL], (0.091, 0.081, 0.071, 0.061))

        weights = t.load_weights()
        self.assertEqual(len(weights), 5)

        self.assertEqual(weights[SNPModelDB.WDBQF.RSID], (u'rs1', u'rs2', u'rs3', u'rs4', u'rs5', u'rs6', u'rs1'))
        self.assertEqual(weights[SNPModelDB.WDBQF.GENE], (u'A', u'A', u'A', u'B', u'B', u'C', u'D'))
        self.assertEqual(weights[SNPModelDB.WDBQF.WEIGHT], (0.2, 0.1, 0.05, 0.4, 0.3, 0.5, 0.6))
        self.assertEqual(weights[SNPModelDB.WDBQF.REF_ALLELE], (u'C', u'A', u'G', u'T', u'C', u'T', u'T'))
        self.assertEqual(weights[SNPModelDB.WDBQF.EFF_ALLELE], (u'T', u'G', u'A', u'C', u'T', u'C', u'C'))

    def test_snps_in_db(self):
        expected = {"rs245915", "rs245913", "rs245909", "rs245906", "rs10486599", "rs144012121", "rs117887801", "rs542000", "rs544632",
                    "rs498475", "rs849327", "rs849336", "rs849335", "rs1513272", "rs849135", "rs849134", "rs860262", "rs849133", "rs1635852",
                    "rs864745", "rs112751321", "rs144273091", "rs117462481", "rs149305679", "rs643036", "rs1937888", "rs17155745", "rs62626328"}
        actual = SNPModelDB.snps_in_db("tests/_td/dbs/test_2.db")
        self.assertEqual(actual, expected)

if __name__ == '__main__':
    unittest.main()
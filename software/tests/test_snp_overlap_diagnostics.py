import unittest

from metax.misc import SnpOverlapDiagnostics


class TestCheckSnpOverlapFromCounts(unittest.TestCase):
    def test_good_overlap_no_message(self):
        model_snps = {"rs1", "rs2", "rs3"}
        pct, msg = SnpOverlapDiagnostics.check_snp_overlap_from_counts(model_snps, matched_count=3)
        self.assertEqual(pct, 100.0)
        self.assertIsNone(msg)

    def test_zero_matches_produces_message(self):
        model_snps = {"rs1", "rs2", "rs3", "rs4", "rs5"}
        pct, msg = SnpOverlapDiagnostics.check_snp_overlap_from_counts(model_snps, matched_count=0)
        self.assertEqual(pct, 0.0)
        self.assertIsNotNone(msg)
        self.assertIn("0.00%", msg)
        self.assertIn("--variant_mapping", msg)
        self.assertIn("--liftover", msg)

    def test_empty_model_snps_no_message(self):
        pct, msg = SnpOverlapDiagnostics.check_snp_overlap_from_counts(set(), matched_count=0)
        self.assertEqual(pct, 0.0)
        self.assertIsNone(msg)

    def test_min_pct_threshold(self):
        model_snps = {str(i) for i in range(100)}
        pct, msg = SnpOverlapDiagnostics.check_snp_overlap_from_counts(model_snps, matched_count=5, min_pct=10.0)
        self.assertEqual(pct, 5.0)
        self.assertIsNotNone(msg)


class TestIdFormats(unittest.TestCase):
    def test_classify(self):
        S = SnpOverlapDiagnostics
        self.assertEqual(S.classify_snp_id("rs9442371"), S.RSID)
        self.assertEqual(S.classify_snp_id("chr1_1083182_C_T_b38"), S.VARID)
        self.assertEqual(S.classify_snp_id("10_94140_G_A_b38"), S.VARID)
        self.assertEqual(S.classify_snp_id("chr1_2685931_CATGTGACAGCC_C_b38"), S.VARID)
        self.assertEqual(S.classify_snp_id("1:13380:C:G"), S.CHR_POS)
        self.assertEqual(S.classify_snp_id("chr1:13380"), S.CHR_POS)
        self.assertEqual(S.classify_snp_id("affx-12345"), S.UNRECOGNIZED)

    def test_dominant_format_tolerates_a_minority(self):
        snps = {"rs{}".format(i) for i in range(34)} | {"chr1_{}_A_G_b38".format(i) for i in range(3)}
        self.assertEqual(SnpOverlapDiagnostics.dominant_id_format(snps), SnpOverlapDiagnostics.RSID)

    def test_disagreement(self):
        S = SnpOverlapDiagnostics
        self.assertTrue(S.formats_disagree(S.RSID, S.VARID))
        self.assertFalse(S.formats_disagree(S.VARID, S.VARID))
        self.assertFalse(S.formats_disagree(S.UNRECOGNIZED, S.RSID))


if __name__ == '__main__':
    unittest.main()

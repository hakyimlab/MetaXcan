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


if __name__ == '__main__':
    unittest.main()

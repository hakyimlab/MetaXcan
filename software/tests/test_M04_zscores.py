import unittest

from . import SampleData
from metax import PredictionModel
from metax import MatrixManager
D = MatrixManager.GENE_SNP_COVARIANCE_DEFINITION

from metax.metaxcan import Utilities


def _prediction_model():
    e = SampleData.dataframe_from_extra(SampleData.sample_extra_1())
    w = SampleData.dataframe_from_weights(SampleData.sample_weights_1())
    return PredictionModel.Model(w, e)


def _covariance():
    s = SampleData.dataframe_from_covariance(SampleData.sample_covariance_s_1())
    return MatrixManager.MatrixManager(s, D)


def _context_good_overlap():
    gwas = SampleData.dataframe_from_gwas(SampleData.sample_gwas_data_4())
    model = _prediction_model()
    return Utilities._build_context(model, _covariance(), gwas)


def _context_mismatched():
    gwas = SampleData.dataframe_from_gwas(SampleData.sample_gwas_data_4())
    gwas = gwas.copy()
    # replace rsids with varID-like ids that share no overlap with the model's rsids,
    # simulating e.g. a MASHR (varID-keyed) model matched against an rsID gwas.
    gwas["snp"] = ["chr1_{}_A_G_b38".format(i) for i in range(len(gwas))]
    model = _prediction_model()
    return Utilities._build_context(model, _covariance(), gwas)


class TestCheckSnpOverlap(unittest.TestCase):
    def test_good_overlap_no_message(self):
        model_snps = {"rs1", "rs2", "rs3"}
        gwas_snps = {"rs1", "rs2", "rs3", "rs4"}
        pct, msg = Utilities.check_snp_overlap(model_snps, gwas_snps)
        self.assertEqual(pct, 100.0)
        self.assertIsNone(msg)

    def test_mismatched_overlap_produces_message(self):
        model_snps = {"rs1", "rs2", "rs3", "rs4", "rs5"}
        gwas_snps = {"chr1_1_A_G_b38", "chr1_2_A_G_b38"}
        pct, msg = Utilities.check_snp_overlap(model_snps, gwas_snps)
        self.assertEqual(pct, 0.0)
        self.assertIsNotNone(msg)
        self.assertIn("0.00%", msg)
        self.assertIn("--model_db_snp_key varID", msg)
        self.assertIn("--keep_non_rsid", msg)

    def test_no_message_on_real_context_with_good_overlap(self):
        c = _context_good_overlap()
        pct, msg = Utilities.check_snp_overlap(c.get_model_snps(), c.get_gwas_snps())
        self.assertAlmostEqual(pct, 15.0 / 18.0 * 100.0)
        self.assertIsNone(msg)

    def test_message_on_real_context_with_mismatched_ids(self):
        c = _context_mismatched()
        pct, msg = Utilities.check_snp_overlap(c.get_model_snps(), c.get_gwas_snps())
        self.assertEqual(pct, 0.0)
        self.assertIsNotNone(msg)
        self.assertIn("0.00%", msg)


if __name__ == '__main__':
    unittest.main()

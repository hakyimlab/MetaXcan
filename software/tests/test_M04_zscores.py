import copy
import logging
import os
import unittest

from . import SampleData
from metax import PredictionModel
from metax import MatrixManager
D = MatrixManager.GENE_SNP_COVARIANCE_DEFINITION

from metax.metaxcan import Utilities

import M03_betas
import M04_zscores


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

    def test_partial_genome_coverage_is_not_a_mismatch(self):
        # a gwas restricted to one chromosome overlaps a genome-wide model by a
        # few percent, which is expected rather than broken
        model_snps = {"rs{}".format(i) for i in range(1000)}
        gwas_snps = {"rs{}".format(i) for i in range(30)}
        pct, msg = Utilities.check_snp_overlap(model_snps, gwas_snps)
        self.assertEqual(pct, 3.0)
        self.assertIsNone(msg)

    def test_format_mismatch_warns_above_min_pct(self):
        # the few ids a varID-keyed gwas does match in an rsid-keyed model are
        # the model's own varID fallbacks, which put overlap over min_pct
        model_snps = {"rs{}".format(i) for i in range(34)} | {"chr1_{}_A_G_b38".format(i) for i in range(3)}
        gwas_snps = {"chr1_{}_A_G_b38".format(i) for i in range(100)}
        pct, msg = Utilities.check_snp_overlap(model_snps, gwas_snps)
        self.assertGreater(pct, 1.0)
        self.assertIsNotNone(msg)
        self.assertIn("different", msg)
        self.assertIn("--model_db_snp_key varID", msg)

    def test_same_format_zero_overlap_warns(self):
        # e.g. hg19 positions against an hg38 model: same id format, no match
        model_snps = {"chr1_{}_A_G_b38".format(i) for i in range(100)}
        gwas_snps = {"chr1_{}_A_G_b37".format(i) for i in range(100)}
        pct, msg = Utilities.check_snp_overlap(model_snps, gwas_snps)
        self.assertEqual(pct, 0.0)
        self.assertIsNotNone(msg)
        self.assertIn("genome build", msg)


class TestCheckCovarianceOverlap(unittest.TestCase):
    def test_matching_keys_no_message(self):
        snps = {"chr1_{}_A_G_b38".format(i) for i in range(10)}
        pct, msg = Utilities.check_covariance_overlap(snps, snps)
        self.assertEqual(pct, 100.0)
        self.assertIsNone(msg)

    def test_rsid_model_against_varid_covariance_warns(self):
        model_snps = {"rs{}".format(i) for i in range(10)}
        covariance_snps = {"chr1_{}_A_G_b38".format(i) for i in range(10)}
        pct, msg = Utilities.check_covariance_overlap(model_snps, covariance_snps)
        self.assertEqual(pct, 0.0)
        self.assertIsNotNone(msg)
        self.assertIn("--model_db_snp_key varID", msg)

    def test_streamed_covariance_is_skipped(self):
        pct, msg = Utilities.check_covariance_overlap({"rs1"}, None)
        self.assertIsNone(pct)
        self.assertIsNone(msg)


_QGT = os.path.join("tests", "_td", "qgt_chr1_subset")
_GWAS_MESSAGE_MARKERS = ("of the model's SNPs were found in the GWAS data",
                         "No GWAS variants at all survived loading")
_COVARIANCE_MESSAGE_MARKER = "of the model's SNPs were found in the covariance data"


class _SPrediXcanArgs(object):
    """Stands in for SPrediXcan.py's argparse namespace, which M03 and M04 share."""
    def __init__(self, snp_column, model_db_snp_key=None, keep_non_rsid=False):
        self.gwas_folder = os.path.join(_QGT, "gwas")
        self.gwas_file = None
        self.gwas_file_pattern = None
        self.snp_column = snp_column
        self.effect_allele_column = "effect_allele"
        self.non_effect_allele_column = "non_effect_allele"
        self.zscore_column = "zscore"
        self.chromosome_column = None
        self.position_column = None
        self.freq_column = None
        self.beta_column = None
        self.beta_sign_column = None
        self.se_column = None
        self.or_column = None
        self.pvalue_column = None
        self.separator = None
        self.skip_until_header = None
        self.handle_empty_columns = False
        self.input_pvalue_fix = 1e-30
        self.keep_non_rsid = keep_non_rsid
        self.snp_map_file = None
        self.output_folder = None
        self.output = None

        self.model_db_path = os.path.join(_QGT, "mashr_Whole_Blood_chr1_subset.db")
        self.model_db_snp_key = model_db_snp_key
        self.covariance = os.path.join(_QGT, "mashr_Whole_Blood_chr1_subset.txt.gz")
        self.stream_covariance = False
        self.single_snp_model = False
        self.output_file = None
        self.additional_output = False
        self.remove_ens_version = False
        self.overwrite = False
        self.MAX_R = None
        self.gwas_h2 = None
        self.gwas_N = None
        self.verbosity = logging.CRITICAL
        self.throw = True


class TestCheckSnpOverlapRealData(unittest.TestCase):
    """Exercises the diagnostic against a real MASHR model and GWAS.

    Fixtures are a 15-gene chr1 slice of GTEx v8 MASHR Whole Blood (PredictDB)
    and the matching variants of the imputed CARDIoGRAM_C4D_CAD GWAS, both from
    the QGT course data. MASHR models are varID-keyed and the GWAS carries both
    an rsid and a varID column, which is what makes the id-format mismatch that
    the diagnostic exists to catch reproducible here.
    """

    def _run_spredixcan(self, args):
        with self.assertLogs(level=logging.WARNING) as captured:
            M03_args = copy.copy(args)
            M03_args.output_folder = None
            M03_args.output = None
            gwas = M03_betas.run(M03_args)
            results = M04_zscores.run(args, gwas)
        messages = [r.getMessage() for r in captured.records]
        gwas_warnings = [m for m in messages if any(k in m for k in _GWAS_MESSAGE_MARKERS)]
        covariance_warnings = [m for m in messages if _COVARIANCE_MESSAGE_MARKER in m]
        return results, gwas_warnings, covariance_warnings

    def test_mashr_varid_invocation_does_not_warn(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id", model_db_snp_key="varID", keep_non_rsid=True)
        results, gwas_warnings, covariance_warnings = self._run_spredixcan(args)
        self.assertEqual(gwas_warnings, [])
        self.assertEqual(covariance_warnings, [])
        self.assertEqual(len(results), 15)
        self.assertEqual(results.zscore.notnull().sum(), 15)

    def test_varid_gwas_against_rsid_keyed_model_warns(self):
        # the classic MASHR mistake: overlap is 8%, high enough to clear a
        # percentage floor, so the id format is what gives it away
        args = _SPrediXcanArgs(snp_column="panel_variant_id", keep_non_rsid=True)
        results, gwas_warnings, _ = self._run_spredixcan(args)
        self.assertEqual(len(gwas_warnings), 1)
        self.assertIn("8.11%", gwas_warnings[0])
        self.assertIn("--model_db_snp_key varID", gwas_warnings[0])
        self.assertLess(len(results), 15)

    def test_varid_gwas_dropped_as_non_rsid_warns(self):
        # every gwas row is filtered out for not looking like an rsID, leaving
        # nothing to compare formats against
        args = _SPrediXcanArgs(snp_column="panel_variant_id", model_db_snp_key="varID")
        _, gwas_warnings, _ = self._run_spredixcan(args)
        self.assertEqual(len(gwas_warnings), 1)
        self.assertIn("No GWAS variants at all survived loading", gwas_warnings[0])
        self.assertIn("--keep_non_rsid", gwas_warnings[0])

    def test_rsid_keyed_model_against_varid_covariance_warns(self):
        # gwas matching is fine here, but the MASHR covariance is varID keyed,
        # so every gene is computed from zero snps and every zscore is NA
        args = _SPrediXcanArgs(snp_column="variant_id")
        results, gwas_warnings, covariance_warnings = self._run_spredixcan(args)
        self.assertEqual(gwas_warnings, [])
        self.assertEqual(len(covariance_warnings), 1)
        self.assertIn("--model_db_snp_key varID", covariance_warnings[0])
        self.assertEqual(results.n_snps_used.sum(), 0)


if __name__ == '__main__':
    unittest.main()

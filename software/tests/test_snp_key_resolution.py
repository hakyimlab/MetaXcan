import os
import tempfile
import unittest

import pandas

from metax.misc import SnpKeyResolution
from metax.misc import SnpOverlapDiagnostics

from .test_M04_zscores import _SPrediXcanArgs, _QGT

import SPrediXcan

RSID = SnpOverlapDiagnostics.RSID
VARID = SnpOverlapDiagnostics.VARID
UNRECOGNIZED = SnpOverlapDiagnostics.UNRECOGNIZED

_MASHR_COLUMNS = {"rsid": RSID, "varID": VARID}


class TestChooseModelSnpKey(unittest.TestCase):
    def test_switches_to_the_column_matching_the_gwas(self):
        # the classic MASHR mistake: a varID-keyed gwas against the default rsid key
        key, msg = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "rsid", VARID, VARID)
        self.assertEqual(key, "varID")
        self.assertIn("varID", msg)

    def test_switches_when_no_covariance_is_known(self):
        key, _ = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "rsid", VARID, None)
        self.assertEqual(key, "varID")

    def test_leaves_a_working_key_alone(self):
        key, msg = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "varID", VARID, VARID)
        self.assertIsNone(key)
        self.assertIsNone(msg)

    def test_gwas_and_covariance_disagreeing_cannot_be_fixed_by_a_key(self):
        # rsid gwas, varID covariance: matching one empties the other
        key, msg = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "rsid", RSID, VARID)
        self.assertIsNone(key)
        self.assertIn("no single column", msg)

    def test_does_not_switch_when_that_would_break_the_covariance(self):
        key, msg = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "rsid", VARID, RSID)
        self.assertIsNone(key)
        self.assertIn("no single column", msg)

    def test_says_nothing_about_ids_it_cannot_classify(self):
        key, msg = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "rsid", UNRECOGNIZED, VARID)
        self.assertIsNone(key)
        self.assertIsNone(msg)

    def test_no_alternative_column_to_offer(self):
        # an older db with only an rsid column, against a varID gwas
        key, msg = SnpKeyResolution.choose_model_snp_key({"rsid": RSID}, "rsid", VARID, None)
        self.assertIsNone(key)
        self.assertIsNone(msg)

    def test_unknown_current_key_is_left_to_fail_downstream(self):
        key, msg = SnpKeyResolution.choose_model_snp_key(_MASHR_COLUMNS, "not_a_column", VARID, VARID)
        self.assertIsNone(key)
        self.assertIsNone(msg)


class TestShouldKeepNonRsid(unittest.TestCase):
    def test_rsid_gwas_needs_no_change(self):
        self.assertFalse(SnpKeyResolution.should_keep_non_rsid(RSID))

    def test_varid_gwas_would_be_dropped_wholesale(self):
        self.assertTrue(SnpKeyResolution.should_keep_non_rsid(VARID))

    def test_unrecognized_is_left_alone(self):
        self.assertFalse(SnpKeyResolution.should_keep_non_rsid(UNRECOGNIZED))


class TestSampling(unittest.TestCase):
    """Reads the real MASHR fixtures, since the point is to tell formats apart in the wild."""

    def setUp(self):
        self.model_db = os.path.join(_QGT, "mashr_Whole_Blood_chr1_subset.db")
        self.covariance = os.path.join(_QGT, "mashr_Whole_Blood_chr1_subset.txt.gz")
        self.gwas = os.path.join(_QGT, "gwas", "cardiogram_chr1_subset.txt.gz")

    def test_finds_both_id_columns_of_a_mashr_db(self):
        columns = SnpKeyResolution.model_db_id_columns(self.model_db)
        self.assertEqual(sorted(columns), ["rsid", "varID"])

    def test_classifies_each_model_column(self):
        varids = SnpKeyResolution.sample_model_ids(self.model_db, "varID")
        self.assertEqual(SnpOverlapDiagnostics.dominant_id_format(varids), VARID)
        rsids = SnpKeyResolution.sample_model_ids(self.model_db, "rsid")
        # a mashr db's rsid column falls back to the varID string where no rsID
        # exists, so the format is only dominant, never uniform
        self.assertEqual(SnpOverlapDiagnostics.dominant_id_format(rsids), RSID)

    def test_classifies_each_gwas_column(self):
        ids = SnpKeyResolution.sample_gwas_ids(self.gwas, "panel_variant_id")
        self.assertEqual(SnpOverlapDiagnostics.dominant_id_format(ids), VARID)
        ids = SnpKeyResolution.sample_gwas_ids(self.gwas, "variant_id")
        self.assertEqual(SnpOverlapDiagnostics.dominant_id_format(ids), RSID)

    def test_classifies_the_covariance(self):
        ids = SnpKeyResolution.sample_covariance_ids(self.covariance)
        self.assertEqual(SnpOverlapDiagnostics.dominant_id_format(ids), VARID)

    def test_missing_or_unreadable_input_is_not_an_error(self):
        self.assertEqual(SnpKeyResolution.model_db_id_columns("nope.db"), [])
        self.assertIsNone(SnpKeyResolution.sample_covariance_ids("nope.txt.gz"))
        self.assertIsNone(SnpKeyResolution.sample_gwas_ids(self.gwas, "no_such_column"))


class TestResolveSnpKeyArguments(unittest.TestCase):
    def test_classic_mashr_mistake_is_corrected(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id")
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertEqual(args.model_db_snp_key, "varID")
        self.assertTrue(args.keep_non_rsid)
        self.assertEqual(len(messages), 2)

    def test_correct_invocation_is_left_untouched(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id", model_db_snp_key="varID", keep_non_rsid=True)
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertEqual(args.model_db_snp_key, "varID")
        self.assertTrue(args.keep_non_rsid)
        self.assertEqual(messages, [])

    def test_non_rsid_gwas_gets_the_loader_filter_lifted(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id", model_db_snp_key="varID")
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertEqual(args.model_db_snp_key, "varID")
        self.assertTrue(args.keep_non_rsid)
        self.assertEqual(len(messages), 1)
        self.assertIn("--keep_non_rsid", messages[0])

    def test_rsid_gwas_against_varid_covariance_is_reported_not_papered_over(self):
        args = _SPrediXcanArgs(snp_column="variant_id")
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertIsNone(args.model_db_snp_key)
        self.assertFalse(args.keep_non_rsid)
        self.assertEqual(len(messages), 1)
        self.assertIn("no single column", messages[0])

    def test_an_explicit_key_is_never_overridden(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id", model_db_snp_key="rsid", keep_non_rsid=True)
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertEqual(args.model_db_snp_key, "rsid")
        self.assertEqual(len(messages), 1)
        self.assertIn("as requested", messages[0])

    def test_special_header_handling_opts_out(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id")
        args.skip_until_header = "variant_id"
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertIsNone(args.model_db_snp_key)
        self.assertEqual(messages, [])


class TestSPrediXcanEndToEnd(unittest.TestCase):
    """The payoff: the classic mistake now produces the same results as the correct invocation."""

    def _run(self, **kwargs):
        args = _SPrediXcanArgs(**kwargs)
        d = tempfile.mkdtemp()
        args.output_file = os.path.join(d, "results.csv")
        SPrediXcan.run(args)
        return pandas.read_csv(args.output_file)

    def test_mashr_model_with_varid_gwas_and_no_flags(self):
        results = self._run(snp_column="panel_variant_id")
        self.assertEqual(len(results), 15)
        self.assertEqual(results.zscore.notnull().sum(), 15)

    def test_does_not_leak_the_chosen_key_back_to_the_caller(self):
        # MetaMany loops over model dbs with a single args; a key resolved for
        # one db need not even exist in the next
        args = _SPrediXcanArgs(snp_column="panel_variant_id")
        args.output_file = os.path.join(tempfile.mkdtemp(), "results.csv")
        SPrediXcan.run(args)
        self.assertIsNone(args.model_db_snp_key)
        self.assertFalse(args.keep_non_rsid)

    def test_matches_the_explicitly_correct_invocation(self):
        auto = self._run(snp_column="panel_variant_id")
        explicit = self._run(snp_column="panel_variant_id", model_db_snp_key="varID", keep_non_rsid=True)
        self.assertTrue(auto.equals(explicit))


if __name__ == '__main__':
    unittest.main()

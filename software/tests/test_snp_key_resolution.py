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


class TestChooseKeys(unittest.TestCase):
    def test_switches_to_the_column_matching_the_gwas(self):
        # the classic MASHR mistake: a varID-keyed gwas against the default rsid key
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "rsid", VARID, VARID)
        self.assertEqual(model_key, "varID")
        self.assertIsNone(gwas_key)
        self.assertIn("varID", msg)

    def test_switches_when_no_covariance_is_known(self):
        model_key, gwas_key, _ = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "rsid", VARID, None)
        self.assertEqual(model_key, "varID")
        self.assertIsNone(gwas_key)

    def test_leaves_a_working_key_alone(self):
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "varID", VARID, VARID)
        self.assertIsNone(model_key)
        self.assertIsNone(gwas_key)
        self.assertIsNone(msg)

    def test_rsid_gwas_against_varid_covariance_translates_the_gwas(self):
        # the covariance decides the key, and the gwas is translated to it
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "rsid", RSID, VARID)
        self.assertEqual(model_key, "varID")
        self.assertEqual(gwas_key, "rsid")
        self.assertIn("Translating", msg)

    def test_varid_gwas_against_rsid_covariance_translates_the_gwas(self):
        model_key, gwas_key, _ = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "rsid", VARID, RSID)
        self.assertIsNone(model_key)
        self.assertEqual(gwas_key, "varID")

    def test_covariance_belonging_to_another_model_is_reported(self):
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys({"rsid": RSID}, "rsid", RSID, VARID)
        self.assertIsNone(model_key)
        self.assertIsNone(gwas_key)
        self.assertIn("not the format of any id column", msg)

    def test_says_nothing_about_ids_it_cannot_classify(self):
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "rsid", UNRECOGNIZED, VARID)
        self.assertIsNone(model_key)
        self.assertIsNone(gwas_key)
        self.assertIsNone(msg)

    def test_no_column_in_the_gwas_format(self):
        # an older db with only an rsid column, against a varID gwas
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys({"rsid": RSID}, "rsid", VARID, None)
        self.assertIsNone(model_key)
        self.assertIsNone(gwas_key)
        self.assertIsNone(msg)

    def test_unknown_current_key_is_left_to_fail_downstream(self):
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys(_MASHR_COLUMNS, "not_a_column", VARID, VARID)
        self.assertIsNone(model_key)
        self.assertIsNone(gwas_key)
        self.assertIsNone(msg)

    def test_an_explicit_key_is_honored_and_the_gwas_translated_to_it(self):
        model_key, gwas_key, msg = SnpKeyResolution.choose_keys(
            _MASHR_COLUMNS, "rsid", VARID, VARID, allow_key_switch=False)
        self.assertIsNone(model_key)
        self.assertEqual(gwas_key, "varID")
        self.assertIn("as requested", msg)


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
        self.assertEqual(SnpKeyResolution.model_db_id_columns(os.path.join(tempfile.mkdtemp(), "nope.db")), [])
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

    def test_rsid_gwas_against_varid_covariance_keys_on_the_covariance(self):
        args = _SPrediXcanArgs(snp_column="variant_id")
        messages = SnpKeyResolution.resolve_snp_key_arguments(args)
        self.assertEqual(args.model_db_snp_key, "varID")
        self.assertEqual(args.gwas_snp_key, "rsid")
        self.assertFalse(args.keep_non_rsid)
        self.assertEqual(len(messages), 1)
        self.assertIn("Translating", messages[0])

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

    def test_hand_mapped_ids_opt_out(self):
        args = _SPrediXcanArgs(snp_column="panel_variant_id")
        args.snp_map_file = "some_map.txt"
        self.assertEqual(SnpKeyResolution.resolve_snp_key_arguments(args), [])
        args = _SPrediXcanArgs(snp_column="variant_id")
        args.gwas_snp_key = "rsid"
        self.assertEqual(SnpKeyResolution.resolve_snp_key_arguments(args), [])


class TestModelIdTranslation(unittest.TestCase):
    def setUp(self):
        self.model_db = os.path.join(_QGT, "mashr_Whole_Blood_chr1_subset.db")

    def test_maps_the_db_rsids_onto_its_varids(self):
        t = SnpKeyResolution.model_id_translation(self.model_db, "rsid", "varID")
        self.assertEqual(t["rs374313793"], "chr1_1689221_G_A_b38")
        # the rsid column falls back to the varID string where no rsID exists
        self.assertEqual(t["chr1_1704673_G_A_b38"], "chr1_1704673_G_A_b38")
        self.assertEqual(set(t.values()), set(SnpKeyResolution.sample_model_ids(self.model_db, "varID")))

    def test_translating_a_column_to_itself_is_nothing(self):
        self.assertEqual(SnpKeyResolution.model_id_translation(self.model_db, "rsid", "rsid"), {})

    def test_unreadable_db_is_not_an_error(self):
        missing = os.path.join(tempfile.mkdtemp(), "nope.db")
        self.assertEqual(SnpKeyResolution.model_id_translation(missing, "rsid", "varID"), {})
        self.assertEqual(SnpKeyResolution.sample_model_ids(missing, "rsid"), [])
        # sqlite3.connect would have created one
        self.assertFalse(os.path.exists(missing))


class TestSPrediXcanEndToEnd(unittest.TestCase):
    """The payoff: the classic mistake now produces the same results as the correct invocation."""

    def _run(self, gwas_folder=None, **kwargs):
        args = _SPrediXcanArgs(**kwargs)
        if gwas_folder:
            args.gwas_folder = gwas_folder
        d = tempfile.mkdtemp()
        args.output_file = os.path.join(d, "results.csv")
        SPrediXcan.run(args)
        return pandas.read_csv(args.output_file)

    def _allele_flipped_gwas_folder(self):
        """The fixture GWAS with every variant reported against the other allele.

        Swapping the effect and non-effect alleles and negating the zscore
        describes the same association, so it has to produce the same result.
        This is the path where a mistake would be silent and wrong rather than
        loud and empty, and it has to hold with an id translation in the way.
        """
        gwas = pandas.read_table(os.path.join(_QGT, "gwas", "cardiogram_chr1_subset.txt.gz"))
        flipped = gwas.assign(effect_allele=gwas.non_effect_allele,
                              non_effect_allele=gwas.effect_allele,
                              zscore=-gwas.zscore)
        folder = tempfile.mkdtemp()
        # na_rep because the loader splits on whitespace, and an empty field
        # would shift every column after it on the ten rows that have one
        flipped.to_csv(os.path.join(folder, "flipped.txt.gz"), sep="\t", index=False,
                       na_rep="NA", compression="gzip")
        return folder

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

    def test_rsid_gwas_against_a_varid_keyed_covariance(self):
        # the case a column choice can't fix: matching the gwas on rsid leaves
        # the varID-keyed covariance empty, so the gwas is translated instead
        results = self._run(snp_column="variant_id")
        self.assertEqual(len(results), 15)
        self.assertEqual(results.zscore.notnull().sum(), 15)

    def test_translated_gwas_agrees_with_the_untranslated_one(self):
        # the same variants named two ways in the same gwas file. Three of the
        # 37 can't be translated -- the db has no rsID for them, its rsid column
        # holds the varID string instead -- so those genes lose a snp. Every
        # gene that kept all of its snps has to come out bit for bit the same:
        # a translation that flipped an allele or mangled a value would show up
        # here as a sign change, not as a missing snp.
        # results come out sorted by pvalue, which the two runs need not agree on
        translated = self._run(snp_column="variant_id").set_index("gene").sort_index()
        direct = self._run(snp_column="panel_variant_id", model_db_snp_key="varID",
                           keep_non_rsid=True).set_index("gene").sort_index()

        same_snps = translated.n_snps_used == direct.n_snps_used
        self.assertEqual(same_snps.sum(), 12)
        self.assertTrue(translated[same_snps].equals(direct[same_snps]))
        # and the rest differ only by having had fewer snps to work with
        self.assertTrue((translated.n_snps_used <= direct.n_snps_used).all())
        self.assertEqual(translated.n_snps_used.sum(), direct.n_snps_used.sum() - 3)

    def test_allele_flip_survives_the_translation(self):
        folder = self._allele_flipped_gwas_folder()
        plain = self._run(snp_column="variant_id").set_index("gene").sort_index()
        flipped = self._run(snp_column="variant_id", gwas_folder=folder).set_index("gene").sort_index()
        self.assertEqual(flipped.n_snps_used.sum(), plain.n_snps_used.sum())
        self.assertTrue(flipped.equals(plain))

    def test_allele_flip_survives_the_key_switch(self):
        folder = self._allele_flipped_gwas_folder()
        plain = self._run(snp_column="panel_variant_id").set_index("gene").sort_index()
        flipped = self._run(snp_column="panel_variant_id", gwas_folder=folder).set_index("gene").sort_index()
        self.assertEqual(flipped.n_snps_used.sum(), plain.n_snps_used.sum())
        self.assertTrue(flipped.equals(plain))


if __name__ == '__main__':
    unittest.main()

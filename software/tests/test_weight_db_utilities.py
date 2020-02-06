import unittest
import sys

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.WeightDBUtilities as WeightDBUtilities


class TestWeightDBUtilities(unittest.TestCase):
    def testGeneEntry(self):
        entry = WeightDBUtilities.GeneEntry("a", "b", "c", "d", "e", "f")
        self.assertEqual(entry.gene, "a")
        self.assertEqual(entry.gene_name, "b")
        self.assertEqual(entry.n_snps, "c")
        self.assertEqual(entry.pred_perf_R2, "d")
        self.assertEqual(entry.pred_perf_pval, "e")
        self.assertEqual(entry.pred_perf_qval, "f")

    def testWeightDBEntry(self):
        entry = WeightDBUtilities.WeightDBEntry("a", "b", "c", "d", "e")
        self.assertEqual(entry.rsid, "a")
        self.assertEqual(entry.gene, "b")
        self.assertEqual(entry.weight, "c")
        self.assertEqual(entry.ref_allele, "d")
        self.assertEqual(entry.eff_allele, "e")

    def testWeightDBInvalidPath(self):
        weight_db = WeightDBUtilities.WeightDB("tests/kk.db")

        with self.assertRaises(RuntimeError):
            weight_db.openDBIfNecessary()

    def testWeightDB(self):
        #test setup
        class DummyCallback():
            def __init__(self):
                self.entries = []

            def __call__(self, weight, extra):
                self.entries.append((weight, extra))

        expected_weights = expected_weights_results()
        expected_extra = expected_extra_results()

        weight_db = WeightDBUtilities.WeightDB("tests/_td/test.db")

        #load gene data
        extra = weight_db.loadExtraColumnData("A")
        self.assertExtra(extra, [expected_extra[0]])

        extra = weight_db.loadExtraColumnData("B")
        self.assertExtra(extra, [expected_extra[1]])

        extra = weight_db.loadExtraColumnData("C")
        self.assertExtra(extra, [expected_extra[2]])

        extra = weight_db.loadExtraColumnData("D")
        self.assertExtra(extra, [expected_extra[3]])

        extra = weight_db.loadExtraColumnData()
        self.assertExtra(extra, expected_extra)

        #load db
        callback = DummyCallback()
        weights = weight_db.loadFromDB(callback, "A")
        self.assertWeights(weights, [expected_weights[0], expected_weights[1], expected_weights[2]])
        self.assertEqual(len(callback.entries),3)
        callback_weights = [e[0] for e in callback.entries]
        self.assertEqual(callback_weights, weights)

        callback = DummyCallback()
        weights = weight_db.loadFromDB(callback, "B")
        self.assertWeights(weights, [expected_weights[3], expected_weights[4]])
        self.assertEqual(len(callback.entries),2)
        callback_weights = [e[0] for e in callback.entries]
        self.assertEqual(callback_weights, weights)

        callback = DummyCallback()
        weights = weight_db.loadFromDB(callback, "C")
        self.assertWeights(weights, [expected_weights[5]])
        self.assertEqual(len(callback.entries),1)
        callback_weights = [e[0] for e in callback.entries]
        self.assertEqual(callback_weights, weights)

        callback = DummyCallback()
        weights = weight_db.loadFromDB(callback, "D")
        self.assertWeights(weights, [expected_weights[6]])
        self.assertEqual(len(callback.entries),1)
        callback_weights = [e[0] for e in callback.entries]
        self.assertEqual(callback_weights, weights)

        callback = DummyCallback()
        weights = weight_db.loadFromDB(callback)
        self.assertWeights(weights, expected_weights)
        self.assertEqual(len(callback.entries),7)
        callback_weights = [e[0] for e in callback.entries]
        self.assertEqual(callback_weights, weights)

        #gene names
        gene_names = weight_db.loadGeneNamesFromDB()
        self.assertEqual(gene_names, ["A", "B", "C", "D"])

    def testWeightDBEntryLogic(self):
        weight_db_entry_logic = WeightDBUtilities.WeightDBEntryLogic("tests/_td/test.db")

        expected_weights = expected_weights_results()
        expected_extra = expected_extra_results()

        self.assertEqual(len(weight_db_entry_logic.weights_by_gene), len(expected_extra))
        self.assertEqual(len(weight_db_entry_logic.gene_data_for_gene), len(expected_extra))

        for e in expected_extra:
            self.assertTrue(e.gene in weight_db_entry_logic.weights_by_gene)
            self.assertTrue(e.gene in weight_db_entry_logic.gene_data_for_gene)

            actual_gene_data = weight_db_entry_logic.gene_data_for_gene[e.gene]
            self.assertExtra([actual_gene_data], [e])

            actual_weights = [w for k,w in weight_db_entry_logic.weights_by_gene[e.gene].items()]
            e_w = [w for w in expected_weights if w.gene == e.gene]
            self.assertWeights(actual_weights, e_w)

        self.assertEqual(len(weight_db_entry_logic.genes_for_an_rsid), 6)
        for rsid, genes in weight_db_entry_logic.genes_for_an_rsid.items():
            expected = [w.gene for w in expected_weights if w.rsid == rsid]
            self.assertEqual(expected, genes)

    def assertWeights(self, weights, expected):
        self.assertEqual(len(weights), len(expected))
        for i,actual in enumerate(weights):
            e = expected[i]
            self.assertEqual(actual.rsid, e.rsid)
            self.assertEqual(actual.gene, e.gene)
            self.assertEqual(actual.weight, e.weight)
            self.assertEqual(actual.ref_allele, e.ref_allele)
            self.assertEqual(actual.eff_allele, e.eff_allele)

    def assertExtra(self, extra, expected):
        self.assertEqual(len(extra), len(expected))

        for i,actual in enumerate(extra):
            e = expected[i]
            self.assertEqual(e.gene, actual.gene)
            self.assertEqual(e.gene_name, actual.gene_name)
            self.assertEqual(e.n_snps, actual.n_snps)
            self.assertEqual(e.pred_perf_R2, actual.pred_perf_R2)
            self.assertEqual(e.pred_perf_pval, actual.pred_perf_pval)
            self.assertEqual(e.pred_perf_qval, actual.pred_perf_qval)

def expected_weights_results():
    class DummyWeight(object):
        pass

    weights = []
    expected_data = [
        ["rs1", "A", 0.2, "C", "T"],
        ["rs2", "A", 0.1, "A", "G"],
        ["rs3", "A", 0.05, "G", "A"],
        ["rs4", "B", 0.4, "T", "C"],
        ["rs5", "B", 0.3, "C", "T"],
        ["rs6", "C", 0.5, "T", "C"],
        ["rs1", "D", 0.6, "T", "C"]
    ]

    for e in expected_data:
        w = DummyWeight()
        w.rsid, w.gene, w.weight, w.ref_allele, w.eff_allele = e[0], e[1], e[2], e[3], e[4]
        weights.append(w)

    return weights

def expected_extra_results():
    class DummyExtra(object):
        pass

    extra = []

    expected_data = [
        ["A", "gene1", 3, 0.9, 0.09, 0.091],
        ["B", "gene2", 2, 0.8, 0.08, 0.081],
        ["C", "gene3", 1, 0.7, 0.07, 0.071],
        ["D", "gene4", 1, 0.6, 0.06, 0.061]
    ]

    for e in expected_data:
        entry = DummyExtra()
        entry.gene, entry.gene_name, entry.n_snps, entry.pred_perf_R2, entry.pred_perf_pval, entry.pred_perf_qval = e[0], e[1], e[2], e[3], e[4], e[5]
        extra.append(entry)

    return extra
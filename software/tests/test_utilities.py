import unittest
import sys
import re

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.Utilities as Utilities
import metax.Exceptions as Exceptions

class TestUtilities(unittest.TestCase):
    def testHapName(self):
        hap_name = Utilities.hapName("a")
        self.assertEqual(hap_name, "a.hap.gz")

    def testLegendName(self):
        legend_name = Utilities.legendName("a")
        self.assertEqual(legend_name, "a.legend.gz")

    def testDosageName(self):
        dosage_name = Utilities.dosageName("a")
        self.assertEqual(dosage_name, "a.dosage.gz")

    def testDosageNamesFromFolder(self):
        names = Utilities.dosageNamesFromFolder("tests/_td/dosage_set_1")
        self.assertEqual(names, [])

    def testLegendNamesFromFolder(self):
        names = Utilities.legendNamesFromFolder("tests/_td/dosage_set_1")
        self.assertEqual(names, ["set_chr1"])

    def testHapNamesFromFolder(self):
        names = Utilities.hapNamesFromFolder("tests/_td/dosage_set_1")
        self.assertEqual(names, ["set_chr1"])

    def testNamesWithPatternFromFolders(self):
        names = Utilities.namesWithPatternFromFolder("tests/_td/dosage_set_1/", ".sample")
        self.assertEqual(names, ["set"])

    def testContentsWithPatternsFromFolders(self):
        contents = Utilities.contentsWithPatternsFromFolder("tests/_td/dosage_set_1", ["sample", "Fail"])
        contents = {c for c in contents}
        self.assertEqual(contents, set([]))

        contents = Utilities.contentsWithPatternsFromFolder("tests/_td/dosage_set_1", ["set", "sample"])
        contents = {c for c in contents}
        self.assertEqual(contents, {"set.sample"})

    def testContentsWithRegexpFromFolder(self):
        contents = Utilities.contentsWithRegexpFromFolder("tests/_td/dosage_set_1", re.compile(".*sample"))
        self.assertEqual(contents, ["set.sample"])

    def testSamplesInputPath(self):
        path = Utilities.samplesInputPath("tests/_td/dosage_set_1")
        self.assertEqual(path, "tests/_td/dosage_set_1/set.sample")

    def testCheckSubdirectorySanity(self):
        b = Utilities.checkSubdirectorySanity("tests", "tests")
        self.assertFalse(b)

        b = Utilities.checkSubdirectorySanity("tests", "tests/_td")
        self.assertTrue(b)

        b = Utilities.checkSubdirectorySanity("tests/_td", "tests")
        self.assertFalse(b)

    def testFileIterator(self):
        class DummyCallback():
            def __init__(self):
                self.lines = []

            def __call__(self, i, line):
                self.lines.append((i, line.strip()))

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set.sample", header="a")
        with self.assertRaises(Exceptions.MalformedInputFile):
            f.iterate(c)

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set.sample")
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, "ID POP GROUP SEX"),
             (1, "ID1 K HERO male"),
             (2, "ID2 K HERO female"),
             (3, "DI5 K HERO male"),
             (4, "ID3 K HERO female"),
             (5,"B1 L T female")])

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set.sample", "")
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, "ID1 K HERO male"),
             (1, "ID2 K HERO female"),
             (2, "DI5 K HERO male"),
             (3, "ID3 K HERO female"),
             (4,"B1 L T female")])

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set.sample", "ID POP GROUP SEX")
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, "ID1 K HERO male"),
             (1, "ID2 K HERO female"),
             (2, "DI5 K HERO male"),
             (3, "ID3 K HERO female"),
             (4,"B1 L T female")])

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set.sample", "DI5 K", ignore_until_header=True)
        f.iterate(c)
        self.assertEqual(c.lines,
             [(0, "ID3 K HERO female"),
             (1,"B1 L T female")])

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set_chr1.legend.gz", header="a", compressed=True)
        with self.assertRaises(Exceptions.MalformedInputFile):
            f.iterate(c)

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set_chr1.legend.gz", compressed=True)
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, "id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL"),
             (1, "1:10177:A:AC 10177 A AC Biallelic_INDEL 0.490922844175492 0.360230547550432 0.336309523809524 0.405566600397614 0.494887525562372 0.425319488817891"),
             (2, "rs1:1:A:T 10505 A T Biallelic_SNP 0 0 0 0 0 0"),
             (3, "1:12:C:G 10506 C G Biallelic_SNP 0 0 0 0 0 0"),
             (4, "rs2:2:G:A 10511 G A Biallelic_SNP 0 0 0 0 0 0"),
             (5, "rs3:3:C:T 10511 G A Biallelic_SNP 0 0 0 0 0 0"),
             (6, "rs4:4:C:T 10511 G A Biallelic_SNP 0 0 0 0 0 0")]
        )

        c = DummyCallback()
        f = Utilities.FileIterator("tests/_td/dosage_set_1/set_chr1.legend.gz", header="id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL", compressed=True)
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, "1:10177:A:AC 10177 A AC Biallelic_INDEL 0.490922844175492 0.360230547550432 0.336309523809524 0.405566600397614 0.494887525562372 0.425319488817891"),
             (1, "rs1:1:A:T 10505 A T Biallelic_SNP 0 0 0 0 0 0"),
             (2, "1:12:C:G 10506 C G Biallelic_SNP 0 0 0 0 0 0"),
             (3, "rs2:2:G:A 10511 G A Biallelic_SNP 0 0 0 0 0 0"),
             (4, "rs3:3:C:T 10511 G A Biallelic_SNP 0 0 0 0 0 0"),
             (5, "rs4:4:C:T 10511 G A Biallelic_SNP 0 0 0 0 0 0")]
        )

    def testCSVFileIterator(self):
        class DummyCallback():
            def __init__(self):
                self.lines = []

            def __call__(self, i, row):
                self.lines.append((i, row))

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set.sample", header="a")
        with self.assertRaises(Exceptions.MalformedInputFile):
            f.iterate(c)

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set.sample")
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, ["ID", "POP", "GROUP", "SEX"]),
             (1, ["ID1", "K", "HERO", "male"]),
             (2, ["ID2", "K", "HERO", "female"]),
             (3, ["DI5", "K", "HERO", "male"]),
             (4, ["ID3", "K", "HERO", "female"]),
             (5, ["B1", "L", "T", "female"])]
        )

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set.sample", "")
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, ["ID1", "K", "HERO", "male"]),
             (1, ["ID2", "K", "HERO", "female"]),
             (2, ["DI5", "K", "HERO", "male"]),
             (3, ["ID3", "K", "HERO", "female"]),
             (4, ["B1", "L", "T", "female"])]
        )

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set.sample", "DI5 K", ignore_until_header=True)
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, ["ID3", "K", "HERO", "female"]),
             (1, ["B1", "L", "T", "female"])]
        )

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set.sample", header="ID POP GROUP SEX")
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, ["ID1", "K", "HERO", "male"]),
             (1, ["ID2", "K", "HERO", "female"]),
             (2, ["DI5", "K", "HERO", "male"]),
             (3, ["ID3", "K", "HERO", "female"]),
             (4, ["B1", "L", "T", "female"])]
        )


        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set_chr1.legend.gz", header="a", compressed=True)
        with self.assertRaises(Exceptions.MalformedInputFile):
            f.iterate(c)

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set_chr1.legend.gz", compressed=True)
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, ["id", "position", "a0", "a1", "TYPE", "AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]),
             (1, ["1:10177:A:AC", "10177", "A", "AC", "Biallelic_INDEL", "0.490922844175492", "0.360230547550432", "0.336309523809524", "0.405566600397614", "0.494887525562372", "0.425319488817891"]),
             (2, ["rs1:1:A:T", "10505", "A", "T", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (3, ["1:12:C:G", "10506", "C", "G", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (4, ["rs2:2:G:A", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (5, ["rs3:3:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (6, ["rs4:4:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"])]
        )

        c = DummyCallback()
        f = Utilities.CSVFileIterator("tests/_td/dosage_set_1/set_chr1.legend.gz", header="id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL", compressed=True)
        f.iterate(c)
        self.assertEqual(c.lines,
            [(0, ["1:10177:A:AC", "10177", "A", "AC", "Biallelic_INDEL", "0.490922844175492", "0.360230547550432", "0.336309523809524", "0.405566600397614", "0.494887525562372", "0.425319488817891"]),
             (1, ["rs1:1:A:T", "10505", "A", "T", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (2, ["1:12:C:G", "10506", "C", "G", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (3, ["rs2:2:G:A", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (4, ["rs3:3:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
             (5, ["rs4:4:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"])]
        )

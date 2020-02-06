import unittest
import sys
import shutil
import os
import io
import gzip
import re

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.ThousandGenomesUtilities as ThousandGenomesUtilities

def buildDummyPeople():
    people_info = [
        "ID1 K HERO male",
        "ID2 K HERO female",
        "DI5 K HERO male",
        "ID3 K HERO female",
        "B1 L T female"
    ]

    people = []

    class Dummy(object):
        pass

    for info in people_info:
        p = Dummy()
        components = info.split()
        p.id, p.population, p.group, p.sex = components[0], components[1], components[2], components[3]
        people.append(p)

    return people

def setupData(root):
    if os.path.exists(root):
        shutil.rmtree(root)
    shutil.copytree("tests/_td/dosage_set_1", os.path.join(root,"dosage_set_1"))
    shutil.copy("tests/_td/snp.txt.gz", root)

def cleanUpData(root):
    shutil.rmtree(root)

class TestThousandGenomesUtilities(unittest.TestCase):
    def testLegendLoader(self):
        class DummyCalback:
            def __init__(self):
                self.lines = []

            def __call__(self, index, row):
                self.lines.append((index, row))

        callback = DummyCalback()
        loader = ThousandGenomesUtilities.LEGENDLoader("tests/_td/dosage_set_1", "set_chr1")
        loader.iterateOverFileLegends(callback)
        self.assertEqual(6, len(callback.lines))
        expected = [
            (0, ["1:10177:A:AC", "10177", "A", "AC", "Biallelic_INDEL", "0.490922844175492", "0.360230547550432", "0.336309523809524", "0.405566600397614", "0.494887525562372", "0.425319488817891"]),
            (1, ["rs1:1:A:T", "10505", "A", "T", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (2, ["1:12:C:G", "10506", "C", "G", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (3, ["rs2:2:G:A", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (4, ["rs3:3:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"]),
            (5, ["rs4:4:C:T", "10511", "G", "A", "Biallelic_SNP", "0", "0", "0", "0", "0", "0"])
        ]

        for i, line in enumerate(callback.lines):
            e = expected[i]
            self.assertEqual(e, line)

    def testImputeLoader(self):
        class DummyCalback:
            def __init__(self):
                self.lines = []

            def __call__(self, hap_line, legend_line):
                self.lines.append((legend_line.strip(), hap_line.strip()))

        callback = DummyCalback()
        loader = ThousandGenomesUtilities.IMPUTELoader("tests/_td/dosage_set_1", "set_chr1")
        loader.iterateOverFile(callback)
        self.assertEqual(6, len(callback.lines))
        expected = [
            ("1:10177:A:AC 10177 A AC Biallelic_INDEL 0.490922844175492 0.360230547550432 0.336309523809524 0.405566600397614 0.494887525562372 0.425319488817891", "0 0 0 0 0 0 0 0 0 0"),
            ("rs1:1:A:T 10505 A T Biallelic_SNP 0 0 0 0 0 0", "0 1 0 1 1 1 1 1 0 0"),
            ("1:12:C:G 10506 C G Biallelic_SNP 0 0 0 0 0 0", "1 1 1 1 1 1 1 1 1 1"),
            ("rs2:2:G:A 10511 G A Biallelic_SNP 0 0 0 0 0 0", "0 0 0 1 0 0 0 0 0 0"),
            ("rs3:3:C:T 10511 G A Biallelic_SNP 0 0 0 0 0 0", "0 0 0 1 1 0 1 0 0 0"),
            ("rs4:4:C:T 10511 G A Biallelic_SNP 0 0 0 0 0 0", "1 1 1 1 1 1 1 1 1 1")
        ]

        for i, line in enumerate(callback.lines):
            e = expected[i]
            self.assertEqual(e, line)

    def testFilteredDosageFileBuilderNoOutput(self):
        builder = ThousandGenomesUtilities.IMPUTEFilteredDosageFileBuilder(
            base_path="_test/dosage_set_1", name="set_chr1", output_pattern="_test/result/set_chr1", chromosome_name="chr1")

        #nothing should be picked up
        setupData("_test")
        os.mkdir("_test/result")
        builder.buildPrediXcan()
        self.assertPredixcanOutput("_test/result/set_chr1.dosage.gz", [])
        cleanUpData("_test")

    def testFilteredDosageFileBuilder(self):
        all_people = buildDummyPeople()
        builder = ThousandGenomesUtilities.IMPUTEFilteredDosageFileBuilder(
            base_path="_test/dosage_set_1", name="set_chr1", output_pattern="_test/result/set_chr1", chromosome_name="chr1")
        builder.all_people = all_people
        builder.selected_people_by_id = {p.id:p for p in all_people}

        #A problematic snp should be discarded
        builder.snp_dict = {"rs1":True}
        setupData("_test")
        os.mkdir("_test/result")
        builder.buildPrediXcan()
        self.assertPredixcanOutput("_test/result/set_chr1.dosage.gz",
        [])
        cleanUpData("_test")

        #Everything in a single snp should be picked up
        builder.snp_dict = {"rs2":True}
        setupData("_test")
        os.mkdir("_test/result")
        builder.buildPrediXcan()
        self.assertPredixcanOutput("_test/result/set_chr1.dosage.gz",
            ["chr1 rs2 2 G A 0.1 0 1 0 0 0"])
        cleanUpData("_test")

        #Filtering snp info
        builder.snp_dict = {"rs3":True, "rs4":True}
        builder.selected_people_by_id = {all_people[0].id:all_people[0], all_people[1].id:all_people[1] }
        setupData("_test")
        os.mkdir("_test/result")
        builder.buildPrediXcan()
        self.assertPredixcanOutput("_test/result/set_chr1.dosage.gz",
        ["chr1 rs3 3 G A 0.25 0 1",
         "chr1 rs4 4 G A 1.0 2 2"])
        cleanUpData("_test")


    def assertPredixcanOutput(self, file, expected):
        lines = []
        with gzip.open(file) as f:
            for l in f:
                lines.append(l.decode().strip())
        self.assertEqual(len(expected), len(lines))

        for i,l in enumerate(lines):
            e = expected[i]
            self.assertEqual(l, e)
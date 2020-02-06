#!/usr/bin/env python


import sys
import pandas
import numpy
import numpy.testing
# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest

from metax.genotype import GTExGenotype
from metax import Utilities

def compare_data_frames(dataframe, dataframe_2, gtex_ids):
    numpy.testing.assert_array_equal(dataframe.rsid, dataframe_2.RS_ID_dbSNP142_CHG37p13)
    numpy.testing.assert_array_equal(dataframe.chromosome, dataframe_2.Chr)
    numpy.testing.assert_array_equal(dataframe.position, dataframe_2.Pos)
    numpy.testing.assert_array_equal(dataframe.ref_allele, dataframe_2.Ref_b37)
    numpy.testing.assert_array_equal(dataframe.alt_allele, dataframe_2.Alt)
    for id in gtex_ids:
        numpy.testing.assert_array_almost_equal(dataframe[id], dataframe_2[id], err_msg="gtex geno failed for %s " % (id,))

class TestGTExGenotype(unittest.TestCase):
    "Each test for save also tests the corresponding load function"
    def test_gtex_geno_lines_generator(self):
        data = []
        for i, line in enumerate(GTExGenotype.gtex_geno_lines("tests/_td/genotype/gtex_like.txt.gz", "tests/_td/genotype/gtex_snp.txt.gz")):
            data.append(line)

        header = GTExGenotype.gtex_geno_header("tests/_td/genotype/gtex_like.txt.gz")
        gtex_ids = header[1:]
        header = ["rsid", "chromosome", "position", "ref_allele", "alt_allele", "frequency"]+gtex_ids
        dataframe = Utilities.to_dataframe(data, header, to_numeric="ignore")


        gtex_snp = pandas.read_table("tests/_td/genotype/gtex_snp.txt.gz")
        dataframe_2 = pandas.read_table("tests/_td/genotype/gtex_like.txt.gz")
        dataframe_2 = pandas.merge(dataframe_2,gtex_snp, left_on="Id", right_on="VariantID")

        compare_data_frames(dataframe, dataframe_2, gtex_ids)


    def test_gtex_geno_by_chromosome(self):
        header = GTExGenotype.gtex_geno_header("tests/_td/genotype/gtex_like.txt.gz")
        gtex_ids = header[1:]

        gtex_snp = pandas.read_table("tests/_td/genotype/gtex_snp.txt.gz")
        dataframe_2 = pandas.read_table("tests/_td/genotype/gtex_like.txt.gz")
        dataframe_2 = pandas.merge(dataframe_2,gtex_snp, left_on="Id", right_on="VariantID")

        def torture_dataframe(dataframe_2, chromosome):
            d = dataframe_2.loc[dataframe_2.Chr == chromosome] if chromosome is not None else dataframe_2
            d = pandas.DataFrame(d)
            rsids = list(d.RS_ID_dbSNP142_CHG37p13)
            used = set()
            for rsid in rsids:
                if rsid in used:
                    i = d[d.RS_ID_dbSNP142_CHG37p13 == rsid].index.tolist()[0]
                    d = d.drop(i)
                used.add(rsid)
            d["number"] = list(range(0, len(d)))
            d = d.set_index("number")
            return d

        def torture_dosage(metadata, dosage, gtex_ids):
            d = [dosage[x] for x in metadata.rsid]
            d = Utilities.to_dataframe(d, gtex_ids, to_numeric="ignore")
            d["rsid"] = list(metadata.rsid)
            d = pandas.merge(metadata, d, on="rsid")
            d["number"] = list(range(0, len(d)))
            d = d.set_index("number")
            return d

        m = pandas.DataFrame()
        for chromosome, metadata, dosage in GTExGenotype.gtex_geno_by_chromosome("tests/_td/genotype/gtex_like.txt.gz", "tests/_td/genotype/gtex_snp.txt.gz"):
            m = pandas.concat([m, metadata])

            self.assertEqual(len(list(dosage.keys())), metadata.shape[0])
            for key in list(dosage.keys()):
                self.assertEqual(len(dosage[key]), 100) #this data has 100 samples
            d = torture_dosage(metadata, dosage, gtex_ids)
            d2 = torture_dataframe(dataframe_2, chromosome)
            compare_data_frames(d, d2, gtex_ids)

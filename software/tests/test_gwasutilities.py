#!/usr/bin/env python

import sys
import gzip
import os
import numpy.random
import metax.Exceptions as Exceptions
from tests import silentRm

# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.GWASUtilities as GWASUtilities
import unittest

class GwasEntry(object):
    def __init__(self, chr, pos, rsid):
        self.chr = chr
        self.pos = pos
        self.rsid = rsid
        self.allele1 = numpy.random.choice(list("ACGT"))
        self.allele2 = numpy.random.choice(list(set(list("ACGT")) - set(self.allele1)))
        self.stderr, self.pval, self.freq = numpy.random.rand(3)
        self.oddsr, self.beta, self.betas, self.betaz = numpy.random.rand(4)

    def __str__(self):
        return "\t".join([str(x) for x in [self.chr, self.pos, self.rsid, self.allele1,
                                 self.allele2, self.stderr, self.pval, self.freq,
                                self.oddsr, self.beta, self.betas, self.betaz]])

class TestGWASUtilities(unittest.TestCase):
    def setUp(self):
        self.filename = "__gwasutil_wheader__.txt"
        self.header_line = "Position,Marker,Allele1,Allele2,StdErr,P-value,Freq1,OddsRatio,Beta,BetaSign,BetaZscore".split(",")

        gzfile = gzip.open("%s.gz" % (self.filename), 'w')
        with open(self.filename, 'w') as file:
            print >> file, "\t".join(self.header_line)
            print >> gzfile, "\t".join(self.header_line)

            self.gwas_entries = []
            for i in range(1, 10):
                self.gwas_entries.append(GwasEntry('1', 1000 * i,
                                "rs%d " % (100000 + numpy.random.randn())))
                print >> file, str(self.gwas_entries[-1])
                print >> gzfile, str(self.gwas_entries[-1])

        gzfile.close()


    def tearDown(self):
        os.remove("%s.gz" % (self.filename))
        os.remove(self.filename)

    def testGWASUtilitiesAddColNoErrors(self):
        gwff = GWASUtilities.GWASFileFormat(self.filename, compressed=False)
        gwff.addSNPColumn("Marker")
        gwff.addSEColumn("StdErr")
        gwff.addA1Column("Allele1")
        gwff.addA2Column("Allele2")
        gwff.addPValueColumn("P-value")
        gwff.addFrequencyColumn("Freq1")
        gwff.addORColumn("OddsRatio")
        gwff.addBetaColumn("Beta")
        gwff.addBetaSignColumn("BetaSign")
        gwff.addBetaZScoreColumn("BetaZscore")
        self.assertEqual(gwff.SNP, 1)
        self.assertEqual(gwff.A1, 2)
        self.assertEqual(gwff.A2, 3)
        self.assertEqual(gwff.SE, 4)
        self.assertEqual(gwff.P, 5)
        self.assertEqual(gwff.FRQ, 6)
        self.assertEqual(gwff.OR, 7)
        self.assertEqual(gwff.BETA, 8)
        self.assertEqual(gwff.BETA_SIGN, 9)
        self.assertEqual(gwff.BETA_ZSCORE, 10)
    def testGWASUtilitiesAddColNoErrorsCompressed(self):
        gwff = GWASUtilities.GWASFileFormat("%s.gz" % self.filename, compressed=True)
        gwff.addSNPColumn("Marker")
        gwff.addSEColumn("StdErr")
        gwff.addA1Column("Allele1")
        gwff.addA2Column("Allele2")
        gwff.addPValueColumn("P-value")
        gwff.addFrequencyColumn("Freq1")
        gwff.addORColumn("OddsRatio")
        gwff.addBetaColumn("Beta")
        gwff.addBetaSignColumn("BetaSign")
        gwff.addBetaZScoreColumn("BetaZscore")
        self.assertEqual(gwff.SNP, 1)
        self.assertEqual(gwff.A1, 2)
        self.assertEqual(gwff.A2, 3)
        self.assertEqual(gwff.SE, 4)
        self.assertEqual(gwff.P, 5)
        self.assertEqual(gwff.FRQ, 6)
        self.assertEqual(gwff.OR, 7)
        self.assertEqual(gwff.BETA, 8)
        self.assertEqual(gwff.BETA_SIGN, 9)
        self.assertEqual(gwff.BETA_ZSCORE, 10)

    def testGWASUtilitiesAddColErrorsCompressed(self):
        filename = "%s.gz" % self.filename
        gwff = GWASUtilities.GWASFileFormat(filename, compressed=True)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addSNPColumn("Markerzz")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("SNP", "Markerzz", filename)), fnexcpt.exception.msg)


        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addSEColumn("StdErrzz")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("SE", "StdErrzz", filename)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addA1Column("all1")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("A1", "all1", filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addA2Column("all2")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("A2", "all2", filename)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addFrequencyColumn("F")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("frequency", "F", filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addPValueColumn("P")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("pvalue", "P", filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addORColumn("OR")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("OR", "OR", filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaColumn("Betas")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("beta", "Betas", filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaSignColumn("BetaSigns")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("beta sign", "BetaSigns", filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaZScoreColumn("BetaZcore")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("beta zscore", "BetaZcore", filename)), fnexcpt.exception.msg)


    def testGWASUtilitiesAddColErrors(self):
        gwff = GWASUtilities.GWASFileFormat(self.filename, compressed=False)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addSNPColumn("Markerzz")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("SNP", "Markerzz", self.filename)), fnexcpt.exception.msg)


        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addSEColumn("StdErrzz")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("SE", "StdErrzz", self.filename)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addA1Column("all1")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("A1", "all1", self.filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addA2Column("all2")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("A2", "all2", self.filename)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addFrequencyColumn("F")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("frequency", "F", self.filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addPValueColumn("P")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("pvalue", "P", self.filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addORColumn("OR")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("OR", "OR", self.filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaColumn("Betas")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("beta", "Betas", self.filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaSignColumn("BetaSigns")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("beta sign", "BetaSigns", self.filename)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaZScoreColumn("BetaZcore")
        self.assertEqual(("%s column name -%s- not found in file '%s'" % ("beta zscore", "BetaZcore", self.filename)), fnexcpt.exception.msg)


if __name__ == "__main__":
    unittest.main()

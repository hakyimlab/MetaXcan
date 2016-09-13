#!/usr/bin/env python

import sys
import gzip
import os
import numpy.random
#from tests import silentRm

# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import metax.Exceptions as Exceptions
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

    def to_str(self, sep):
        return sep.join([str(x) for x in [self.chr, self.pos, self.rsid, self.allele1,
                                 self.allele2, self.stderr, self.pval, self.freq,
                                self.oddsr, self.beta, self.betas, self.betaz]])

    HEADER_COMPS = "Position,Marker,Allele1,Allele2,StdErr,P-value,Freq1,OddsRatio,Beta,BetaSign,BetaZscore".split(",")

class TestGWASUtilities(unittest.TestCase):
    def setUp(self):
        self.tab_file_name = "__gwasutil_wheader__.txt"

        with gzip.open("%s.gz" % (self.tab_file_name), 'w') as gzfile:
            with open(self.tab_file_name, 'w') as file:
                header_line = "%s\n" % ("\t".join(GwasEntry.HEADER_COMPS))
                file.write(header_line)
                gzfile.write(header_line)

                self.gwas_entries = []
                for i in range(1, 10):
                    self.gwas_entries.append(GwasEntry('1', 1000 * i,
                                "rs%d " % (100000 + numpy.random.randn())))
                    line = "%s\n" % (self.gwas_entries[-1].to_str("\t"))
                    file.write(line)
                    gzfile.write(line)

        self.comma_file_name = "__gwasutil_wheader_c__.txt"

        with gzip.open("%s.gz" % (self.comma_file_name), 'w') as gzfile:
            with open(self.comma_file_name, 'w') as file:
                header_line = "%s\n" % (",".join(GwasEntry.HEADER_COMPS))
                file.write(header_line)
                gzfile.write(header_line)

                self.gwas_entries = []
                for i in range(1, 10):
                    self.gwas_entries.append(GwasEntry('1', 1000 * i,
                                                   "rs%d " % (100000 + numpy.random.randn())))
                    line = "%s\n" % (self.gwas_entries[-1].to_str(","))
                    file.write(line)
                    gzfile.write(line)

    def tearDown(self):
        os.remove("%s.gz" % (self.tab_file_name))
        os.remove(self.tab_file_name)
        os.remove("%s.gz" % (self.comma_file_name))
        os.remove(self.comma_file_name)

#Expectation utilities
    def expectFileFormat(self, f):
        self.assertEqual(f.SNP, 1)
        self.assertEqual(f.OTHER_ALLELE, 3)
        self.assertEqual(f.EFFECT_ALLELE, 2)
        self.assertEqual(f.SE, 4)
        self.assertEqual(f.P, 5)
        self.assertEqual(f.FRQ, 6)
        self.assertEqual(f.OR, 7)
        self.assertEqual(f.BETA, 8)
        self.assertEqual(f.BETA_SIGN, 9)
        self.assertEqual(f.BETA_ZSCORE, 10)

#test methods for file format itself
    def expectGWASUtilitiesAddColNoErrors(self, file_name, compressed, separator=None):
        gwff = GWASUtilities.GWASFileFormat(file_name, compressed=compressed, separator=separator)
        gwff.addSNPColumn("Marker")
        gwff.addSEColumn("StdErr")
        gwff.addOtherAlleleColumn("Allele2")
        gwff.addEffectAlleleColumn("Allele1")
        gwff.addPValueColumn("P-value")
        gwff.addFrequencyColumn("Freq1")
        gwff.addORColumn("OddsRatio")
        gwff.addBetaColumn("Beta")
        gwff.addBetaSignColumn("BetaSign")
        gwff.addBetaZScoreColumn("BetaZscore")
        self.expectFileFormat(gwff)

    def testGWASUtilitiesAddColNoErrors(self):
        self.expectGWASUtilitiesAddColNoErrors(self.tab_file_name, compressed=False)

    def testGWASUtilitiesAddColNoErrorsCompressed(self):
        file_name = "%s.gz" % self.tab_file_name
        self.expectGWASUtilitiesAddColNoErrors(file_name, compressed=True)

    def testGWASUtilitiesAddColNoErrorsTab(self):
        self.expectGWASUtilitiesAddColNoErrors(self.tab_file_name, compressed=False, separator="\t")

    def testGWASUtilitiesAddColNoErrorsCompressedTab(self):
        file_name = "%s.gz" % self.tab_file_name
        self.expectGWASUtilitiesAddColNoErrors(file_name, compressed=True, separator="\t")

    def testGWASUtilitiesAddColNoErrorsWithComma(self):
        self.expectGWASUtilitiesAddColNoErrors(self.comma_file_name, compressed=False, separator=",")

    def testGWASUtilitiesAddColNoErrorsCompressedWithComma(self):
        file_name = "%s.gz" % self.comma_file_name
        self.expectGWASUtilitiesAddColNoErrors(file_name, compressed=True, separator=",")

    def expectGWASUtilitiesAddColErrors(self, file_name, compressed, separator=None):
        gwff = GWASUtilities.GWASFileFormat(file_name, compressed=compressed, separator=separator)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addSNPColumn("Markerzz")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("SNP", "Markerzz", file_name)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addSEColumn("StdErrzz")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("SE", "StdErrzz", file_name)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addOtherAlleleColumn("all2")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("-other allele-", "all2", file_name)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addEffectAlleleColumn("all1")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("-effect allele-", "all1", file_name)), fnexcpt.exception.msg)

        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addFrequencyColumn("F")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("frequency", "F", file_name)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addPValueColumn("P")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("pvalue", "P", file_name)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addORColumn("OR")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("OR", "OR", file_name)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaColumn("Betas")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("beta", "Betas", file_name)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaSignColumn("BetaSigns")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("beta sign", "BetaSigns", file_name)), fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.ReportableException) as fnexcpt:
            gwff.addBetaZScoreColumn("BetaZcore")
        self.assertEqual(("%s column name -%s- not found in file '%s'. Is the file compressed?" % ("beta zscore", "BetaZcore", file_name)), fnexcpt.exception.msg)

    def testGWASUtilitiesAddColErrorsCompressed(self):
        filename = "%s.gz" % self.tab_file_name
        self.expectGWASUtilitiesAddColErrors(filename, compressed=True)

    def testGWASUtilitiesAddColErrors(self):
        self.expectGWASUtilitiesAddColErrors(self.tab_file_name, compressed=False)

#test methods for file foramt factory
    def expectGWASFileFormatFromArgs(self, file_name, compressed, separator=None):
        class Dummy(object): pass
        args = Dummy()

        #
        #Invalid arguments shouldn't work. Test them one at a time
        with self.assertRaises(AttributeError) as fnexcpt:
            f = GWASUtilities.GWASFileFormat.fileFormatFromArgs(file_name, args)

        args.separator = separator
        with self.assertRaises(AttributeError) as fnexcpt:
            f = GWASUtilities.GWASFileFormat.fileFormatFromArgs(file_name, args)

        args.compressed = compressed
        with self.assertRaises(AttributeError) as fnexcpt:
            f = GWASUtilities.GWASFileFormat.fileFormatFromArgs(file_name, args)

        args.skip_until_header = False
        with self.assertRaises(AttributeError) as fnexcpt:
            f = GWASUtilities.GWASFileFormat.fileFormatFromArgs(file_name, args)

        args.snp_column="Marker"
        args.frequency_column = "Freq1"
        args.effect_allele_column = "Allele1"
        args.other_allele_column = "Allele2"
        args.pvalue_column = "P-value"
        args.or_column = "OddsRatio"
        args.frequency_column = "Freq1"
        args.beta_column = "Beta"
        args.beta_sign_column = "BetaSign"
        args.beta_zscore_column = "BetaZscore"
        args.se_column = "StdErr"

        #"Position,Marker,Allele1,Allele2,StdErr,P-value,Freq1,OddsRatio,Beta,BetaSign,BetaZscore"
        f = GWASUtilities.GWASFileFormat.fileFormatFromArgs(file_name, args)
        self.expectFileFormat(f)

    def testGWASFileFormatFromArgs(self):
        self.expectGWASFileFormatFromArgs(self.tab_file_name, compressed=False)

    def testGWASFileFormatFromArgsCompressed(self):
        file_name = "%s.gz" % self.tab_file_name
        self.expectGWASFileFormatFromArgs(file_name, compressed=True)

    def testGWASFileFormatFromArgsTab(self):
        self.expectGWASFileFormatFromArgs(self.tab_file_name, compressed=False, separator="\t")

    def testGWASFileFormatFromArgsTabCompressed(self):
        file_name = "%s.gz" % self.tab_file_name
        self.expectGWASFileFormatFromArgs(file_name, compressed=True, separator="\t")

if __name__ == "__main__":
    unittest.main()

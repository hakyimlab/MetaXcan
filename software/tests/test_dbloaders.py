#!/usr/bin/env python

import sys
import os
import sqlite3      # Just use direct DBI to build up the database
import metax.Exceptions as Exceptions
import numpy
from . import silentRm

# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest
from metax.deprecated.DBLoaders import DBLoaders

class geneEntry(object):
    def __init__(self, gene, genename, rsq, snp_count):
        self.gene = gene
        self.genename = genename
        self.rsq = float(rsq)
        self.snp_count = int(snp_count)

    def commit(self, dbcrs):
        dbcrs.execute("INSERT INTO genes VALUES (?,?,?,?)",
                      (self.gene, self.genename, self.rsq, self.snp_count,))
# cursor.execute("SELECT rsid1, rsid2, covariance FROM covariances")
class covarEntry(object):
    def __init__(self, rsid1, rsid2, covar):
        self.rsid1 = rsid1
        self.rsid2 = rsid2
        self.covar = covar

    def commit(self, dbcrs):
        dbcrs.execute("INSERT INTO covariances VALUES (?,?,?)",
                      (self.rsid1, self.rsid2, self.covar))


class varianceEntry(object):
    def __init__(self, rsid, var):
        self.rsid = rsid
        self.var = var

    def commit(self, dbcrs):
        dbcrs.execute("insert into variances VAlUES (?,?)",
                      (self.rsid, self.var))

class TestDBLoaderBase(unittest.TestCase):
    def setUp(self):
        self.db_filename = "__dbloader_test__.db"

        # Sqlite doesn't like to create tables if they already exist. This weird filename should be safe enough...
        silentRm(self.db_filename)

        cnx = sqlite3.connect(self.db_filename)
        crs = cnx.cursor()

        self.genes = [
            geneEntry("Gene1", "G1B", "0.2142", "7"),
            geneEntry("Gene2", "GB2", "0.5124", "13"),
            geneEntry("Three", "R3GN", "0.001", "2"),
            geneEntry("Gene4", "ABCB", "0.6", "5")
        ]

        crs.execute("CREATE TABLE genes (gene TEXT, genename TEXT, R2 DOUBLE, 'snp_count' INTEGER)")
        for gene in self.genes:
            gene.commit(crs)
        cnx.commit()


        self.covars = [
            covarEntry("rs1", "rs2", 0.5),
            covarEntry("rs1", "rs1", 0.1),
            covarEntry("rs1", "rs3", 0.3),
            covarEntry("rs1", "rs4", 0.4),
            covarEntry("rs2", "rs2", 0.2),
            covarEntry("rs2", "rs3", 0.7),
            covarEntry("rs2", "rs7", 0.8),
            covarEntry("rs3", "rs8", 0.8),
            covarEntry("rs3", "rs3", 0.9)
        ]
        crs.execute("CREATE TABLE covariances (rsid1 TEXT, rsid2 TEXT, covariance DOUBLE)")
        for covar in self.covars:
            covar.commit(crs)
        cnx.commit()

        self.covariances = numpy.array([0.1, 0.5, 0.3, 0.5, 0.2, 0.7, 0.3, 0.7, 0.9]).reshape(3,3)
        self.cv_valid_keys = ["rs1", "rs2", "rs3"]


        self.variances = [
            varianceEntry("rs1234", 0.5),
            varianceEntry("rs3456", 0.6),
            varianceEntry("rs5432", 0.1),
            varianceEntry("rs64352", 0.2)
        ]
        crs.execute("CREATE TABLE variances (rsid TEXT, var DOUBLE)")
        for var in self.variances:
            var.commit(crs)
        cnx.commit()

    def tearDown(self):
        os.remove(self.db_filename)

class TestDbLoaderBasic(TestDBLoaderBase):
    def testBasics(self):
        rsquareds = DBLoaders.loadKeyedDataSetFromDB(self.db_filename, "genes", "gene", "R2")
        genenames = DBLoaders.loadKeyedDataSetFromDB(self.db_filename, "genes", "gene", "genename")
        snpcount = DBLoaders.loadKeyedDataSetFromDB(self.db_filename, "genes", "gene", "snp_count")
        self.assertEqual(self.db_filename, rsquareds.name)
        idx = 0
        for gene in self.genes:
            self.assertEqual(gene.rsq, rsquareds.values_by_key[gene.gene], "RSquared for %s" % (gene.gene))
            self.assertEqual(gene.genename, genenames.values_by_key[gene.gene], "Genename for %s (%s, %s)" % (gene.gene, gene.genename,genenames.values_by_key[gene.gene]))
            self.assertEqual(gene.snp_count, snpcount.values_by_key[gene.gene], "Snp Count for %s (%s, %s)" % (gene.gene, gene.snp_count, snpcount.values_by_key[gene.gene]))
            idx += 1

    def testVariance(self):
        variances = DBLoaders.loadVariancesFromDB(self.db_filename)
        self.assertEqual(len(self.variances), len(variances.values_by_key))
        for var in self.variances:
            self.assertAlmostEqual(var.var, variances.values_by_key[var.rsid], 0.01)

    def testCovariances(self):
        covariances, valid_keys = DBLoaders.loadCovarianceMatrix(self.db_filename, None)
        self.assertEqual(len(self.cv_valid_keys), len(valid_keys))
        self.assertTrue(numpy.array_equal(self.covariances, covariances))


    def testExceptions(self):
        with self.assertRaises(Exceptions.BadFilename) as fnexcpt:
            badfilename = DBLoaders.loadKeyedDataSetFromDB("bad_filename", "genes", "gene", "R2")
        self.assertEqual("Invalid filename: bad_filename", fnexcpt.exception.msg)
        with self.assertRaises(Exceptions.InvalidDbFormat) as fnexcpt:
            badtablename = DBLoaders.loadKeyedDataSetFromDB(self.db_filename, "tissues", "gene", "R2")
        self.assertEqual("\t'no such table: %s'"% ("tissues"), fnexcpt.exception.msg.split("\n")[1])
        with self.assertRaises(Exceptions.InvalidDbFormat) as fnexcpt:
            badcolone = DBLoaders.loadKeyedDataSetFromDB(self.db_filename, "genes", "tissue", "R2")
        self.assertEqual("\t'no such column: %s'"% ("tissue"), fnexcpt.exception.msg.split("\n")[1])
        with self.assertRaises(Exceptions.InvalidDbFormat) as fnexcpt:
            badcoltwo = DBLoaders.loadKeyedDataSetFromDB(self.db_filename, "genes", "gene", "R1")
        self.assertEqual("\t'no such column: %s'"% ("R1"), fnexcpt.exception.msg.split("\n")[1])




if __name__ == "__main__":
    unittest.main()

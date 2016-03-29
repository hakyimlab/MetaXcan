#!/usr/bin/env python

import sys
import os
import sqlite3      # Just use direct DBI to build up the database

# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")

from metax.DBLoaders import DBLoaders
import unittest

class geneEntry(object):
    def __init__(self, gene, genename, rsq, snp_count):
        self.gene = gene
        self.genename = genename
        self.rsq = float(rsq)
        self.snp_count = int(snp_count)

    def commit(self, dbcrs):
        dbcrs.execute("INSERT INTO genes VALUES (?,?,?,?)",
                      (self.gene, self.genename, self.rsq, self.snp_count,))

class TestDBLoaderBase(unittest.TestCase):
    def setUp(self):
        self.db_filename = "__dbloader_test__.db"

        # Sqlite doesn't like to create tables if they already exist. This weird filename should be safe enough...
        os.system("rm %s "% (self.db_filename))
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


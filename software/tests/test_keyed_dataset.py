#!/usr/bin/env python

import sys
# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")

import os
import gzip
import unittest


from metax.KeyedDataSet import KeyedDataSet
from metax.KeyedDataSet import KeyedDataSetFileUtilities

class TestKeyedDataset(unittest.TestCase):
    def testKeyedDataSet(self):
        keys = ['k1','k2','k3','4','5']
        values = ['v1','v2',100,0.4, 5]

        kds1 = KeyedDataSet("Test1", None, values, keys)
        # Test DataSet construction
        self.assertEqual(kds1.name, "Test1")
        self.assertEqual(kds1.index, None)
        self.assertEqual(len(kds1.data), len(values))
        self.assertEqual(kds1.values_by_key["k1"], 'v1')
        self.assertEqual(kds1.values_by_key["k2"], 'v2')
        self.assertEqual(kds1.values_by_key["k3"], 100)
        self.assertEqual(kds1.values_by_key["4"], 0.4)
        self.assertEqual(kds1.values_by_key["5"], 5)

        kds2 = KeyedDataSet("Test2", 1, keys=keys, data=values)
        self.assertEqual(kds2.name, "Test2")
        self.assertEqual(kds2.index, 1)
        self.assertEqual(kds2.data, values)
        self.assertEqual(len(kds2.data), len(values))

    def testKeyedDataSetWriteContents(self):
        keys = ['k1','k2','k3','4','5']
        values = ['v1','v2',100,0.4, 5]

        kds = KeyedDataSet("Test1", None, values, keys)
        filename = "kds_file.txt"
        with open(filename, 'w') as file:
            KeyedDataSetFileUtilities.writeContents(file, kds, "H1", "H2")
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents,
                         "H1 H2\nk1 v1\nk2 v2\nk3 100\n4 0.4\n5 5")

        # No Header
        with open(filename, 'w') as file:
            KeyedDataSetFileUtilities.writeContents(file, kds, None, None)
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, "k1 v1\nk2 v2\nk3 100\n4 0.4\n5 5")

        os.remove(filename)


    def testKeyedDataSetSave(self):
        keys = ['k1','k2','k3','4','5']
        values = ['v1','v2',100,0.4, 5]

        kds = KeyedDataSet("Test1", None, values, keys)
        filename = "kds1_file.txt"
        KeyedDataSetFileUtilities.saveToFile(filename, kds, "H1", "H2")
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents,
                         "H1 H2\nk1 v1\nk2 v2\nk3 100\n4 0.4\n5 5")

        # No header
        KeyedDataSetFileUtilities.saveToFile(filename, kds, None, None)
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, "k1 v1\nk2 v2\nk3 100\n4 0.4\n5 5")

        os.remove(filename)


    def testKeyedDataSetSaveCompressed(self):
        keys = ['k1','k2','k3','4','5']
        values = ['v1','v2',100,0.4, 5]

        kds = KeyedDataSet("Test1", None, values, keys)
        filename = "kds1_file.txt.gz"
        KeyedDataSetFileUtilities.saveToCompressedFile(filename, kds, "H1", "H2")
        file_contents = gzip.open(filename, 'rb').read().strip()
        self.assertEqual(file_contents,
                         "H1 H2\nk1 v1\nk2 v2\nk3 100\n4 0.4\n5 5")

        # No header
        KeyedDataSetFileUtilities.saveToFile(filename, kds, None, None)
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, "k1 v1\nk2 v2\nk3 100\n4 0.4\n5 5")
        os.remove(filename)

    def testKeyedDataSetsWriteContents(self):
        keys = "k1,k2,k3,k4,k5".split(",")
        values1 = "v1,v2,v3,v4,v5".split(",")
        values2 = "d1,d2,d3,d4,d5".split(",")
        values3 = "j1,j2,j3,j4,j5".split(",")

        kds1 = KeyedDataSet("col1", 1, values1, keys)
        kds2 = KeyedDataSet("col2", 2, values2, keys)
        kds3 = KeyedDataSet("col3", 3, values3, keys)

        filename = "kds_multi.txt"
        with open(filename, "w") as file:
            KeyedDataSetFileUtilities.writeSetsContent(file,
                        [kds1, kds2, kds3], "K")
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, ("K col1 col2 col3\nk1 v1 d1 j1\n" +
                        "k2 v2 d2 j2\nk3 v3 d3 j3\nk4 v4 d4 j4\nk5 v5 d5 j5"))

        with open(filename, "w") as file:
            KeyedDataSetFileUtilities.writeSetsContent(file,
                        [kds1, kds2, kds3], None)
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, ("k1 v1 d1 j1\n" +
                        "k2 v2 d2 j2\nk3 v3 d3 j3\nk4 v4 d4 j4\nk5 v5 d5 j5"))

        os.remove(filename)


    def testKeyedDataSetsSave(self):
        keys = "k1,k2,k3,k4,k5".split(",")
        values1 = "v1,v2,v3,v4,v5".split(",")
        values2 = "d1,d2,d3,d4,d5".split(",")
        values3 = "j1,j2,j3,j4,j5".split(",")

        kds1 = KeyedDataSet("col1", 1, values1, keys)
        kds2 = KeyedDataSet("col2", 2, values2, keys)
        kds3 = KeyedDataSet("col3", 3, values3, keys)

        filename = "kds1_sets_file.txt"
        KeyedDataSetFileUtilities.saveSetsToFile(filename, [kds1, kds2, kds3], "K")
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, ("K col1 col2 col3\nk1 v1 d1 j1\n" +
                        "k2 v2 d2 j2\nk3 v3 d3 j3\nk4 v4 d4 j4\nk5 v5 d5 j5"))
        KeyedDataSetFileUtilities.saveSetsToFile(filename, [kds1, kds2, kds3], None)
        file_contents = open(filename, 'r').read().strip()
        self.assertEqual(file_contents, ("k1 v1 d1 j1\n" +
                        "k2 v2 d2 j2\nk3 v3 d3 j3\nk4 v4 d4 j4\nk5 v5 d5 j5"))

        os.remove(filename)

    def testKeyedDataSetsCompressedSave(self):
        keys = "k1,k2,k3,k4,k5".split(",")
        values1 = "v1,v2,v3,v4,v5".split(",")
        values2 = "d1,d2,d3,d4,d5".split(",")
        values3 = "j1,j2,j3,j4,j5".split(",")

        kds1 = KeyedDataSet("col1", 1, values1, keys)
        kds2 = KeyedDataSet("col2", 2, values2, keys)
        kds3 = KeyedDataSet("col3", 3, values3, keys)

        filename = "kds1_sets_file.txt"
        KeyedDataSetFileUtilities.saveSetsToCompressedFile(filename, [kds1, kds2, kds3], "K")
        file_contents = gzip.open(filename, 'rb').read().strip()
        self.assertEqual(file_contents, ("K col1 col2 col3\nk1 v1 d1 j1\n" +
                        "k2 v2 d2 j2\nk3 v3 d3 j3\nk4 v4 d4 j4\nk5 v5 d5 j5"))
        KeyedDataSetFileUtilities.saveSetsToCompressedFile(filename, [kds1, kds2, kds3], None)
        file_contents = gzip.open(filename, 'rb').read().strip()
        self.assertEqual(file_contents, ("k1 v1 d1 j1\n" +
                        "k2 v2 d2 j2\nk3 v3 d3 j3\nk4 v4 d4 j4\nk5 v5 d5 j5"))

        os.remove(filename)

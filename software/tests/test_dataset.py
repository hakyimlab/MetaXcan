#!/usr/bin/env python

import sys
import os
import numpy

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest
from metax.DataSet import DataSet
from metax.DataSet import DataSetFileUtilities
from metax.DataSet import DataSetCollection

class TestDataset(unittest.TestCase):
    def testDatasetInitialization(self):
        ds1 = DataSet()
        self.assertEqual(0, len(ds1.data))
        self.assertIsNone(ds1.name)
        self.assertIsNone(ds1.index)

        ds2 = DataSet("TEST", 5, [1,2,5,7])
        self.assertEqual(5, ds2.index)
        self.assertEqual("TEST", ds2.name)
        self.assertEqual(4, len(ds2.data))
        self.assertEqual(1, ds2.data[0])
        self.assertEqual(7, ds2.data[-1])

        ds3 = DataSet(name="TEST", index=5, data=[1,2,5,7])
        self.assertEqual(5, ds3.index)
        self.assertEqual("TEST", ds3.name)
        self.assertEqual(4, len(ds3.data))
        self.assertEqual(1, ds3.data[0])
        self.assertEqual(7, ds3.data[-1])

    def testDsfuNoHeader(self, filename="__test_file"):
        with open(filename, "w") as file:
            dataset = numpy.random.randint(10, size=(5,10))
            data = []
            for row in dataset:
                data.append(" ".join([str(x) for x in row]))
                print(data[-1], file=file)
            dataset = data

        ds1 = DataSetFileUtilities.loadFromFile(filename)
        self.assertEqual(filename, ds1.name)
        self.assertEqual(dataset, ds1.data)
        os.remove(filename)


    def testDsfuHeader(self, filename="__test_file"):
        header = "A,B,C,D,E".split(",")
        with open(filename, "w") as file:
            print(" ".join(header), file=file)
            dataset = numpy.random.randint(10, size=(5,10))
            data = []
            for row in dataset:
                data.append(" ".join([str(x) for x in row]))
                print(data[-1], file=file)
            dataset = data
        ds1 = DataSetFileUtilities.loadFromFile(filename, " ".join(header))
        self.assertEqual(filename, ds1.name)
        self.assertEqual(list(dataset), ds1.data)

        os.remove(filename)

if __name__ == "__main__":
    unittest.main()

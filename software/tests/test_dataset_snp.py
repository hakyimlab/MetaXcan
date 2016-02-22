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
from metax.DataSetSNP import DataSetSNP

class TestDatasetSNP(unittest.TestCase):
    def testDatasetInitialization(self):
        ds1 = DataSetSNP()
        self.assertEqual(0, len(ds1.data))
        self.assertIsNone(ds1.name)
        self.assertIsNone(ds1.index)
        self.assertIsNone(ds1.position)
        self.assertIsNone(ds1.ref_allele)
        self.assertIsNone(ds1.eff_allele)

        ds2 = DataSetSNP("TEST", 5, [1,2,5,7], 100, 'A', 'C')
        self.assertEqual(5, ds2.index)
        self.assertEqual("TEST", ds2.name)
        self.assertEqual(4, len(ds2.data))
        self.assertEqual(1, ds2.data[0])
        self.assertEqual(7, ds2.data[-1])
        self.assertEqual(100, ds2.position)
        self.assertEqual('A',ds2.ref_allele)
        self.assertEqual('C', ds2.eff_allele)

        ds3 = DataSetSNP(name="TEST", index=5, data=[1,2,5,7], position=100, ref_allele='A', eff_allele='C')
        self.assertEqual(5, ds3.index)
        self.assertEqual("TEST", ds3.name)
        self.assertEqual(4, len(ds3.data))
        self.assertEqual(1, ds3.data[0])
        self.assertEqual(7, ds3.data[-1])
        self.assertEqual(100, ds3.position)
        self.assertEqual('A',ds3.ref_allele)
        self.assertEqual('C', ds3.eff_allele)

if __name__ == "__main__":
    unittest.main()

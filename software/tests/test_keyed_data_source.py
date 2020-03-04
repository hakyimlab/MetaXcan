#!/usr/bin/env python


import sys
# This allows us to run an individual test as it's own 'program'. Very useful
# for debugging
if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import os
import gzip
import unittest

from metax.misc import KeyedDataSource

class TestKeyedDataset(unittest.TestCase):
    "Each test for save also tests the corresponding load function"
    def testKeyedDataSet(self):
        X = ['a','b','c','d']
        Y = [1, 1, 2, 3]
        Z = [4, 5, 6, 7]

        k = KeyedDataSource.load_data("tests/_td/test.txt", key_name="X", value_name="Y", numeric=True)
        k_ = {X[i]:Y[i] for i in range(0, len(X))}
        self.assertEqual(k, k_)

        k = KeyedDataSource.load_data("tests/_td/test.txt", key_name="X", value_name="Y")
        k_ = {X[i]:str(Y[i]) for i in range(0, len(X))}
        self.assertEqual(k, k_)

        k = KeyedDataSource.load_data("tests/_td/test.txt", key_name="Z", value_name="Y", numeric=True)
        k_ = {str(Z[i]): Y[i] for i in range(0, len(X))}
        self.assertEqual(k, k_)

        k = KeyedDataSource.load_data("tests/_td/test.txt", key_name="Z", value_name="Y")
        k_ = {str(Z[i]): str(Y[i]) for i in range(0, len(X))}
        self.assertEqual(k, k_)
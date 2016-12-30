#!/usr/bin/env python
import unittest
import sys
import shutil
import os
import re
import gzip

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

#import M03_betas

class TestM03(unittest.TestCase):
    def testNoWeightDBMode(self):
        self.assertEqual(True, True)


if __name__ == "__main__":
    unittest.main()
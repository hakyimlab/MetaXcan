#!/usr/bin/env python

import sys
import os

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest




if __name__ == "__main__":
    unittest.main()

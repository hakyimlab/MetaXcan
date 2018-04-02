import numpy
import numpy.testing
import pandas

import unittest
from metax.gwas import GWASSpecialHandling

class TestGWASSpecialHandling(unittest.TestCase):

    def test_sanitize_component(self):
        sanitize = GWASSpecialHandling.sanitize_component
        #Number pass through
        self.assertEqual("a", sanitize("a"))

        #number with comma as decimal separator
        self.assertEqual(".9", sanitize(",9"))
        self.assertEqual("-.9", sanitize("-,9"))
        self.assertEqual("-10.99", sanitize("-10,99"))
        self.assertEqual("10.99", sanitize("10,99"))
        self.assertEqual("+10.99", sanitize("+10,99"))

        self.assertEqual(".9e9", sanitize(",9e9"))
        self.assertEqual(".9e-9", sanitize(",9e-9"))
        self.assertEqual(".9e+9", sanitize(",9e+9"))

        #looks like a number with comma decimal, but is not
        self.assertEqual("a,9", sanitize("a,9"))
        self.assertEqual("-,a9", sanitize("-,a9"))
        self.assertEqual("-10,99e", sanitize("-10,99e"))

        #number with dot as decimal separator
        self.assertEqual(".9", sanitize(".9"))
        self.assertEqual("-.9", sanitize("-.9"))
        self.assertEqual("-10.99", sanitize("-10.99"))
        self.assertEqual("10.99", sanitize("10.99"))
        self.assertEqual("+10.99", sanitize("+10.99"))

        #NA gets translated into Nome
        self.assertIsNone(sanitize("NA"))
        self.assertIsNone(sanitize("."))
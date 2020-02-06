import numpy
import numpy.testing

import unittest
from metax import misc
from metax.misc import FeatureMatrix

from . import SampleData

EXPECTED_A = \
[[ 0.1 ,  0.01,  0.07],
 [ 0.01,  0.02,  0.02],
[ 0.07,  0.02,  0.06]]

EXPECTED_B = \
[[ 3.5  ,  3.84 ,  3.67 ],
[ 3.84 ,  4.44 ,  4.14 ],
[ 3.67 ,  4.14 ,  3.905]]

EXPECTED_D = \
[[ 1.74,  2.44],
[ 2.44,  3.47]]

class TestFeatureMatrix(unittest.TestCase):

    def test_load(self):
        data = SampleData.set_of_feature_sets()
        manager = FeatureMatrix.FeatureMatrixManager(data, standardize=False)

        a = manager.data["a"]
        self.assertEqual(len(a), 3)
        numpy.testing.assert_almost_equal(a["1"], [0.0, 0.1, 0.3, 0.0])
        numpy.testing.assert_almost_equal(a["2"], [0.0, 0.1, 0.0, 0.1])
        numpy.testing.assert_almost_equal(a["3"], [0.0, 0.1, 0.2, 0.1])

        a_m, a_labels = manager.get_feature_product("a")
        numpy.testing.assert_almost_equal(a_m, EXPECTED_A)
        self.assertEqual(a_labels, ["1", "2", "3"])

        a_m, a_labels = manager.get_feature_product("a", center=True)
        numpy.testing.assert_almost_equal(a_m, numpy.cov([a["1"], a["2"], a["3"]]))
        self.assertEqual(a_labels, ["1", "2", "3"])

        b = manager.data["b"]
        self.assertEqual(len(b), 3)
        numpy.testing.assert_almost_equal(b["1"], [1.0, 0.5, 1.2, 0.9])
        numpy.testing.assert_almost_equal(b["2"], [1.0, 1.0, 1.2, 1.0])
        numpy.testing.assert_almost_equal(b["3"], [1.0, 0.75, 1.2, 0.95])

        b_m, b_labels = manager.get_feature_product("b")
        numpy.testing.assert_almost_equal(b_m, EXPECTED_B)
        self.assertEqual(b_labels, ["1", "2", "3"])

        d = manager.data["d"]
        self.assertEqual(len(d), 2)
        numpy.testing.assert_almost_equal(d["2"], [0.5, 0.7, 0.6, 0.8])
        numpy.testing.assert_almost_equal(d["3"], [0.9, 0.9, 0.8, 1.1])

        d_m, d_labels = manager.get_feature_product("d")
        numpy.testing.assert_almost_equal(d_m, EXPECTED_D)
        self.assertEqual(d_labels, ["2", "3"])

        e_m, e_labels = manager.get_feature_product("e")
        numpy.testing.assert_almost_equal(e_m, [[ 0.3726]])
        self.assertEqual(e_labels, ["2"])


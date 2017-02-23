import numpy
import numpy.testing

import unittest
from metax import misc
from metax.misc import FeatureMatrix

import SampleData

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
        manager = FeatureMatrix.FeatureMatrixManager(data)

        a = manager.data["a"]
        self.assertEqual(len(a), 3)
        numpy.testing.assert_almost_equal(a["1"].values, [0.0, 0.1, 0.3, 0.0])
        numpy.testing.assert_almost_equal(a["2"].values, [0.0, 0.1, 0.0, 0.1])
        numpy.testing.assert_almost_equal(a["3"].values, [0.0, 0.1, 0.2, 0.1])

        a_m, a_labels = manager.get_feature_matrix("a")
        numpy.testing.assert_almost_equal(a_m, EXPECTED_A)
        self.assertEqual(a_labels, ["1", "2", "3"])

        b = manager.data["b"]
        self.assertEqual(len(b), 3)
        numpy.testing.assert_almost_equal(b["1"].values, [1.0, 0.5, 1.2, 0.9])
        numpy.testing.assert_almost_equal(b["2"].values, [1.0, 1.0, 1.2, 1.0])
        numpy.testing.assert_almost_equal(b["3"].values, [1.0, 0.75, 1.2, 0.95])

        b_m, b_labels = manager.get_feature_matrix("b")
        numpy.testing.assert_almost_equal(b_m, EXPECTED_B)
        self.assertEqual(b_labels, ["1", "2", "3"])

        d = manager.data["d"]
        self.assertEqual(len(d), 2)
        numpy.testing.assert_almost_equal(d["2"].values, [0.5, 0.7, 0.6, 0.8])
        numpy.testing.assert_almost_equal(d["3"].values, [0.9, 0.9, 0.8, 1.1])

        d_m, d_labels = manager.get_feature_matrix("d")
        numpy.testing.assert_almost_equal(d_m, EXPECTED_D)
        self.assertEqual(d_labels, ["2", "3"])

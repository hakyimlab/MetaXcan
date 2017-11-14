import numpy
import os
import pandas
import numpy.testing

import unittest

from metax.misc import DataFrameStreamer


class TestDataframeStreamer(unittest.TestCase):

    def test_dataframe_streamer(self):
        root = "tests/_td/gtex_like_eqtl"
        kk = os.path.join(root, "data.txt")

        kkk = pandas.DataFrame()
        for d in DataFrameStreamer.data_frame_streamer(kk, "gene_id"):
            gene = d.gene_id.values[0]
            gene = os.path.join(root, gene) + ".txt"
            p = pandas.read_table(gene)
            self.assertTrue(p.equals(d))
            kkk = pandas.concat([kkk, d], ignore_index=True)

        kk = pandas.read_table(kk)
        self.assertTrue(kk.equals(kkk))


if __name__ == "__main__":
    unittest.main()
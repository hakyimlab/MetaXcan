#!/usr/bin/env python
import os
import logging
from metax import Logging, Utilities

import numpy as np
import h5py_cache
import pandas

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists, delete it or move it if you want it generated again")
        return

    Utilities.ensure_requisite_folders(args.output)

    logging.info("Reading input")
    data = pandas.read_table(args.input)

    logging.info("Opening output")
    f = h5py_cache.File(args.output, 'w', chunk_cache_mem_size=int(50 * (1024 ** 2)))

    n_genes = data.shape[1]-2
    n_samples = data.shape[0]
    n_genes_chunk = np.min((n_genes, 10))

    logging.info("Processing expression")
    p = f.create_dataset("pred_expr", shape=(n_genes, n_samples),
                                        chunks=(n_genes_chunk, n_samples),
                                        dtype=np.dtype('float32'), scaleoffset=4, compression='gzip')
    g = f.create_dataset("genes", (n_genes,), dtype="S30")

    for i, gene in enumerate(data.columns.values[2:]):
        p[i, :] = data[gene].to_numpy()
        g[i] = np.string_(gene)

    logging.info("saving samples")
    s = f.create_dataset("samples", (n_samples,), dtype="S25")
    for i in xrange(0, n_samples):
        s[i] = np.string_(data["IID"][i])
    f.close()
    logging.info("Done")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("""Convert expression text file to HDF5.
WARNING: h5py _cache has been deprecated.
Input format should be a tab-separated text file like so:
FID IID GENE1   GENE2 ...
I1  I1  0.1     0.2   ...
I2  I2  0.0     -0.1  ...
....
                                     """)
    parser.add_argument("--input", default="Text file with predicted transcriptomic trait")
    parser.add_argument("--output")
    parser.add_argument("--verbosity", default=logging.INFO, type=int)

    args = parser.parse_args()
    Logging.configureLogging(args.verbosity)
    run(args)
import logging

from MultiPrediXcanAssociation import Context
from .. expression import HDF5Expression

class HDF5Context(Context):
    def __init__(self, args):
        self.args = args
        self.h5 = None

    def __enter__(self):
        logging.info("enter")
        gene_map, h5 = HDF5Expression._structure(self.args.hdf5_expression_folder)
        self.h5 = h5
        self.expression = HDF5Expression.ExpressionManager(gene_map, h5)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logging.info("exit")
        HDF5Expression._close(self.h5)

    def get_genes(self):
        return self.expression.gene_map.keys()

def context_from_args(args):
    logging.info("Preparing context")

    context = HDF5Context(args)

    return context
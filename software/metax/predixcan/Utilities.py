import logging
import pandas

from MultiPrediXcanAssociation import Context
from .. expression import HDF5Expression

class HDF5Context(Context):
    def __init__(self, args):
        self.args = args
        self.h5 = None

    def __enter__(self):
        logging.info("Acquiring HDF5 expression caches")
        gene_map, h5 = HDF5Expression._structure(self.args.hdf5_expression_folder)
        self.h5 = h5
        self.expression = HDF5Expression.ExpressionManager(gene_map, h5)

        logging.info("Accquiring phenotype")
        self.pheno = _pheno_from_file_and_column(self.args.input_phenos_file, self.args.input_phenos_column)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        HDF5Expression._close(self.h5)

    def get_genes(self):
        return self.expression.gene_map.keys()

    def expression_for_gene(self, gene):
        return self.expression.expression_for_gene(gene)

    def get_pheno(self):
        return self.pheno

def _pheno_from_file_and_column(file, column):
    x = pandas.read_table(file)
    p = x[column].values
    return p

def context_from_args(args):
    logging.info("Preparing context")

    context = HDF5Context(args)

    return context
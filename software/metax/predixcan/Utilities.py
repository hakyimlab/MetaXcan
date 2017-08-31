import logging
import pandas
import numpy

from MultiPrediXcanAssociation import Context, MTPMode
from .. expression import HDF5Expression
from .. import Exceptions

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
        if self.args.mode == MTPMode.K_LOGISTIC:
            v = set([str(float(x)) for x in self.pheno])
            if v != {'0.0', '1.0', 'nan'}:
                raise Exceptions.InvalidArguments("Logistic regression was asked but phenotype is not binomial")

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        HDF5Expression._close(self.h5)

    def get_genes(self):
        return self.expression.gene_map.keys()

    def expression_for_gene(self, gene):
        return self.expression.expression_for_gene(gene)

    def get_pheno(self):
        return self.pheno

    def get_mode(self):
        return self.args.mode

def _pheno_from_file_and_column(path, column):
    x = pandas.read_table(path, usecols=[column], sep="\s+")
    p = x[column]
    p.loc[numpy.isclose(p, -999.0, atol=1e-3, rtol=0)] = numpy.nan
    p = p.values
    return p

def _check_args(args):
    if not args.mode in MTPMode.K_MODES:
        raise Exceptions.InvalidArguments("Invalid mode")

def mp_context_from_args(args):
    logging.info("Preparing context")

    _check_args(args)
    context = HDF5Context(args)

    return context
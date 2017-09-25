import logging
import pandas
import numpy

from patsy import dmatrices
import statsmodels.api as sm

from MultiPrediXcanAssociation import Context, MTPMode
from .. expression import HDF5Expression
from .. import Exceptions

class HDF5Context(Context):
    def __init__(self, args):
        self.args = args
        self.h5 = None
        self.mode = None
        self.covariates = None
        self.pheno = None

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

        self.mode = self.args.mode
        if self.args.covariates_file and self.args.covariates:
            self.mode = MTPMode.K_LINEAR
            logging.info("Acquiring covariates")
            self.covariates = _get_covariates(self.args)
            logging.info("Replacing phenotype with residuals")
            self.pheno = _get_residual(self.pheno, self.covariates)

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
        return self.mode

    def get_covariates(self):
        return self.covariates

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

def _get_covariates(args):
    covariates = pandas.read_table(args.covariates_file, sep="\s+")[args.covariates]
    columns = covariates.columns.values
    for c in columns:
        covariates.loc[ numpy.isclose(covariates[c], -999.0, atol=1e-3, rtol=0), c] = numpy.nan
    return covariates

def _get_residual(pheno, covariates):
    e = pandas.DataFrame(covariates)
    e["order"] = xrange(0, e.shape[0])
    e["pheno"] = pheno

    e_ = e.dropna()
    e_ = e_[e_.columns[~(e_ == 0).all()]]

    keys = list(e_.columns.values)
    keys.remove("pheno")
    keys.remove("order")

    y, X = dmatrices("pheno ~ {}".format(" + ".join(keys)), data=e_, return_type="dataframe")
    model = sm.OLS(y,X)
    result = model.fit()
    e_["residual"] = result.resid
    e = e.merge(e_[["order", "residual"]], on="order", how="left")
    return e["residual"]
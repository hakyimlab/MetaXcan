import logging
import pandas
import numpy

from patsy import dmatrices
import statsmodels.api as sm

from .MultiPrediXcanAssociation import Context as _MTPContext, MTPMode
from .PrediXcanAssociation import Context as _PContext, PMode
from metax import Exceptions, Utilities

from .. expression import Expression, PlainTextExpression
try:
    from .. expression import  HDF5Expression
except:
    logging.info("Could not import HDF5 expression")

try:
    import h5py
except:
    logging.info("Could not import h5py module")


from .. import Exceptions

########################################################################################################################
class MTPContext(_MTPContext):
    def __init__(self, args, expression):
        self.args = args
        self.mode = None
        self.covariates = None
        self.pheno = None
        self.pc_filter = _filter_from_args(args)
        self.expression = expression

    def get_genes(self):
        return self.expression.get_genes()

    def expression_for_gene(self, gene):
        return self.expression.expression_for_gene(gene)

    def get_pheno(self):
        return self.pheno

    def get_mode(self):
        return self.mode

    def get_covariates(self):
        return self.covariates

    def get_pc_filter(self):
        return self.pc_filter

    def __enter__(self):
        logging.info("Entering Multi-tissue PrediXcan context")
        self.expression.enter()
        _prepare_phenotype(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logging.info("Exiting Multi-tissue PrediXcan context")
        self.expression.exit()

class DumbMTPContext(_MTPContext):
    def __init__(self, expression, pheno, gene, pc_filter):
        self.expression = expression
        self.pheno = pheno
        self.gene = gene
        self.pc_filter = pc_filter

    def get_genes(self): return [self.gene]
    def expression_for_gene(self, gene): return self.expression
    def get_pheno(self): return self.pheno
    def get_pc_filter(self): return self.pc_filter
    def get_mode(self): return MTPMode.K_LINEAR

def _check_args(args):
    if not args.mode in MTPMode.K_MODES:
        raise Exceptions.InvalidArguments("Invalid mode")

def _expression_manager_from_args(args):
    if args.hdf5_expression_folder:
        logging.info("Preparing expression from HDF5 files")
        expression = HDF5Expression.ExpressionManager(args.hdf5_expression_folder, args.expression_pattern, code_999=args.code_999, standardise=args.standardize_expression)
    elif args.memory_efficient and args.expression_folder:
        logging.info("Loading expression manager from text files (memory efficient)")
        expression = PlainTextExpression.ExpressionManagerMemoryEfficient(args.expression_folder, args.expression_pattern, standardise=args.standardize_expression)
    elif args.expression_folder:
        logging.info("Loading expression manager from text files")
        expression = PlainTextExpression.ExpressionManager(args.expression_folder, args.expression_pattern, standardise=args.standardize_expression)
    else:
        raise RuntimeError("Could not build context from arguments")
    return expression

def mp_context_from_args(args):
    logging.info("Preparing Multi-Tissue PrediXcan context")
    _check_args(args)
    expression = _expression_manager_from_args(args)
    context = MTPContext(args, expression)
    return context

def _filter_eigen_values_from_max(s, ratio):
    s_max = numpy.max(s)
    return [i for i,x in enumerate(s) if x >= s_max*ratio]

def _filter_from_args(args):
    if args.pc_condition_number:
        return lambda x:_filter_eigen_values_from_max(x, 1.0/args.pc_condition_number)
    elif args.pc_eigen_ratio:
        return lambda x:_filter_eigen_values_from_max(x, args.pc_eigen_ratio)
    return None

########################################################################################################################

class PContext(_PContext):
    def __init__(self, args, expression):
        self.args = args
        self.mode = None
        self.covariates = None
        self.pheno = None
        self.expression = expression

    def get_genes(self):
        return self.expression.get_genes()

    def expression_for_gene(self, gene):
        return self.expression.expression_for_gene(gene)

    def get_pheno(self):
        return self.pheno

    def get_mode(self):
        return self.mode

    def get_covariates(self):
        return self.covariates

    def __enter__(self):
        self.expression.enter()
        _prepare_phenotype(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.expression.exit()

class DumbPContext(_PContext):
    def __init__(self, expression, pheno, gene, pc_filter):
        self.expression = expression
        self.pheno = pheno
        self.gene = gene
        self.pc_filter = pc_filter

    def get_genes(self): return [self.gene]
    def expression_for_gene(self, gene): return self.expression
    def get_pheno(self): return self.pheno
    def get_mode(self): return MTPMode.K_LINEAR
    def get_covariates(self): return None

def _check_args_file(args):
    if not args.mode in PMode.K_MODES:
        raise Exceptions.InvalidArguments("Invalid mode")

def expression_from_args(args, prediction_results = None):
    if prediction_results:
        if isinstance(prediction_results, BasicPredictionRepository):
            logging.info("Preparing PrediXcan context from data")
            expression = Expression.ExpressionFromData(prediction_results.genes)
        else:
            raise Exceptions.ReportableException("Invalid prediction results")
    elif args.hdf5_expression_file:
        logging.info("Preparing PrediXcan HDF5 context")
        expression = HDF5Expression.Expression(args.hdf5_expression_file)
    elif args.expression_file:
        logging.info("Preparing PrediXcan context")
        expression = PlainTextExpression.Expression(args.expression_file)
    else:
        raise RuntimeError("Could not build context from arguments")
    return expression

def p_context_from_args(args, prediction_results=None):
    _check_args_file(args)
    expression = expression_from_args(args, prediction_results)
    context = PContext(args, expression)
    return context

########################################################################################################################

def _prepare_phenotype(context):
    logging.info("Accquiring phenotype")
    context.pheno = _pheno_from_file_and_column(context.args.input_phenos_file, context.args.input_phenos_column, context.args.input_phenos_na_values)
    if context.args.mode == MTPMode.K_LOGISTIC:
        try:
            v = set([str(float(x)) for x in context.pheno])
            if not v.issubset({'0.0', '1.0', 'nan'}):
                raise Exceptions.InvalidArguments("Logistic regression was asked but phenotype is not binomial")
        except:
            raise Exceptions.InvalidArguments("Logistic regression: could not parse phenotype")

    context.mode = context.args.mode
    if context.args.covariates_file and context.args.covariates:
        context.mode = MTPMode.K_LINEAR
        logging.info("Acquiring covariates")
        context.covariates = _get_covariates(context.args)
        logging.info("Replacing phenotype with residuals")
        context.pheno = _get_residual(context.pheno, context.covariates)

def _pheno_from_file_and_column(path, column, na_rep=None):
    x = pandas.read_table(path, usecols=[column], sep="\s+", na_values=na_rep)
    p = x[column]
    if not na_rep:
        # Default UKB NA is represented as -999
        p.loc[numpy.isclose(p, -999.0, atol=1e-3, rtol=0)] = numpy.nan
    p = p.values
    return p

def _get_covariates(args):
    covariates = pandas.read_table(args.covariates_file, sep="\s+")[args.covariates]
    columns = covariates.columns.values
    for c in columns:
        covariates.loc[ numpy.isclose(covariates[c], -999.0, atol=1e-3, rtol=0), c] = numpy.nan
    return covariates

def _get_residual(pheno, covariates):
    e = pandas.DataFrame(covariates)
    e["order"] = range(0, e.shape[0])
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

########################################################################################################################

class PredictionRepository:
    def update(self, gene, dosage, weight): raise Exceptions.NotImplemented("gene_repository is not implemented")
    def store_prediction(self): raise Exceptions.NotImplemented("gene_repository is not implemented")
    def summary(self): raise Exceptions.NotImplemented("gene_repository is not implemented")

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

class BasicPredictionRepository(PredictionRepository):
    def __init__(self, samples, extra, output_path):
        self.genes = {}
        self.samples = samples
        self.stats = {}
        self.extra = extra
        self.output_path = output_path

    def update(self, gene, dosage, weight):
        if not gene in self.genes:
            self.stats[gene] = [1, ]
            self.genes[gene] = weight * dosage
        else:
            self.stats[gene][0] += 1
            self.genes[gene] += weight * dosage

    def store_prediction(self):
        logging.info("Saving prediction as a text file")
        d = pandas.DataFrame(self.genes)
        result = pandas.concat([self.samples, d], axis=1, sort=False)
        Utilities.save_dataframe(result, self.output_path)

    def summary(self):
        return summary_report(self.stats, self.extra)
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        pass

class HDF5PredictionRepository(PredictionRepository):
    def __init__(self, samples, extra, path):
        self.path = path
        self.samples = samples
        self.extra = extra
        self.genes_index = None
        self.h5 = None
        self.d_genes = None
        self.d_prediction = None
        self.d_computed = None
        self.d_samples = None
        self.stats = {}
        self.closed = True

    def __enter__(self):
        logging.log(9, "Creating HDF5 prediction repository")
        self.h5 = h5py.File(self.path, 'w')

        n_genes = self.extra.shape[0]
        n_samples = self.samples.shape[0]
        n_genes_chunk = numpy.min((n_genes, 10))

        logging.log(9, "Creating prediction dataset")
        self.d_prediction = self.h5.create_dataset("pred_expr", shape=(n_genes, n_samples), chunks=(n_genes_chunk, n_samples), dtype=numpy.dtype('float32'), scaleoffset=4, compression='gzip')
        self.d_computed = self.h5.create_dataset("computed", (n_genes,), dtype=int)

        logging.log(9, "Creating genes")
        try:
            dt = h5py.string_dtype()
        except:
            logging.info("h5py data initialization failed, last-ditch attempt")
            dt = h5py.special_dtype(vlen=str)

        self.d_genes = self.h5.create_dataset("genes", (n_genes,), dtype=dt)
        self.d_genes[:] = list(map(str, self.extra.gene))
        self.genes_index = {x:i for i,x in enumerate(self.extra.gene)}

        logging.log(9, "Creating samples")
        self.d_samples = self.h5.create_dataset("samples", (n_samples,), dtype=dt)
        self.d_samples[:] = list(map(str,self.samples.IID))

        logging.log(9, "created HDF5 prediction repository")
        self.closed = False
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logging.log(9, "closing HDF5 prediction repository")
        self.h5.file.close()
        self.closed = True
        return self

    def store_prediction(self):
        if not self.closed:
            logging.log(9, "closing HDF5 prediction repository")
            self.h5.file.close()

    def update(self, gene, dosage, weight):
        if not gene in self.stats:
            self.stats[gene] = [1, ]
            self.d_prediction[self.genes_index[gene], :] = weight * dosage
            self.d_computed[self.genes_index[gene]] = 1
        else:
            self.stats[gene][0] += 1
            self.d_prediction[self.genes_index[gene], :] += weight * dosage

    def summary(self):
        return summary_report(self.stats, self.extra)

def summary_report(summary_data, extra):
    s = []
    for k, v in summary_data.items():
        s.append((k, v[0]))
    s = pandas.DataFrame(s, columns=["gene", "n_snps_used"])
    s = extra.merge(s, on="gene", how="left")[['gene', 'gene_name', 'n_snps_in_model', 'n_snps_used', 'pred_perf_r2', 'pred_perf_pval']]
    return s

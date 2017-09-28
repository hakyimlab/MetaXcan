import numpy
import pandas
import logging

from JointAnalysis import Context, ContextMixin
from .. import MatrixManager
from .. import Exceptions
from .. import PredictionModel
from ..misc import GWASAndModels
from ..misc import DataFrameStreamer
from ..misc import KeyedDataSource
from ..genotype import GeneExpressionMatrixManager
from ..metaxcan import MetaXcanResultsManager

class SimpleContext(ContextMixin, Context):
    def __init__(self, metaxcan_results_manager, matrix_manager, genes, cutoff, epsilon):
        self.metaxcan_results_manager = metaxcan_results_manager
        self.matrix_manager = matrix_manager
        self.model_name_index = _index(matrix_manager)

        #Todo: maybe fix?
        self.full_ensemble_id = False
        self.gene_names = self._process_genes(genes)

        # Yeah, a callable object. Go figure.
        self.cutoff = cutoff
        self.epsilon = epsilon

    def get_genes(self):
        matrix_genes = self.matrix_manager.model_labels()
        results_genes = self.metaxcan_results_manager.get_genes()
        genes = {x for x in matrix_genes if x.split(".")[0] in results_genes}
        return genes

    def get_n_genes(self):
        return len(self.get_genes())

    def check(self):
        pass

class ExpressionStreamedContext(ContextMixin, Context):
    """
    A context that works by streaming the snps covariance file, one gene at a time.
    Whenever it yields a gene from `get_genes()`, it will configure itself to that gene,
    but since it holds snp covariance data one gene at a time, it won work if you sort it.
    So, use it as it comes.
    """
    def __init__(self, metaxcan_results_manager, snp_covariance_streamer, model_manager, genes, cutoff, epsilon, full_ensemble_id=False):
        self.metaxcan_results_manager = metaxcan_results_manager
        self.snp_covariance_streamer = snp_covariance_streamer
        self.model_manager = model_manager

        self.full_ensemble_id = full_ensemble_id
        self.gene_names = self._process_genes(genes)

        # Yeah, a callable object. Go figure.
        self.cutoff = cutoff
        self.epsilon = epsilon
        self.matrix_manager = None

    def get_genes(self):
        for d in self.snp_covariance_streamer:
            if not self.full_ensemble_id:
                d.GENE = d.GENE.str.split(".").str.get(0)
            g = d.GENE.values[0]
            self.matrix_manager = GeneExpressionMatrixManager._GeneExpressionMatrixManager(d, self.model_manager)
            yield g

    def get_n_genes(self):
        return len(self.model_manager.get_genes())

    def check(self):
        a = self.model_manager.get_model_labels()
        b = self.metaxcan_results_manager.get_model_labels()
        if len(a.intersection(b)) == 0:
            raise Exceptions.ReportableException("No intersection between model names in MetaXcan Results and Prediction Models. Please verify your input")

def _index(matrix_manager):
    labels = matrix_manager.model_labels()
    labels = {x.split(".")[0]:x for x in labels}
    return labels

def _cutoff(args):
    class CutoffEigenRatio(object):
        def __init__(self, cutoff_ratio):
            self.cutoff_ratio = float(cutoff_ratio)

        def __call__(self, matrix):
            # conceptual shotcut
            if self.cutoff_ratio == 0:
                return 0.0
            w, v = numpy.linalg.eig(matrix)
            w = -numpy.sort(-w)
            cutoff = self.cutoff_ratio * w[0]
            return cutoff

    class CutoffTraceRatio(object):
        def __init__(self, cutoff_ratio):
            self.cutoff_ratio = float(cutoff_ratio)

        def __call__(self, matrix):
            # conceptual shotcut
            if self.cutoff_ratio == 0:
                return 0.0
            trace = numpy.trace(matrix)
            cutoff = self.cutoff_ratio * trace
            return cutoff

    class CutoffThreshold(object):
        def __init__(self, cutoff_threshold):
            self.cutoff_threshold = float(cutoff_threshold)

        def __call__(self, matrix):
            #conceptual shotcut
            if self.cutoff_threshold== 0:
                return 0.0
            eigen = sorted(numpy.linalg.eig(matrix)[0], reverse=True)
            trace = numpy.sum(eigen)
            cumsum = numpy.cumsum(eigen)
            objective = trace*(1-self.cutoff_threshold)
            last = None
            for i in xrange(0, len(cumsum)):
                last = eigen[i]
                if cumsum[i] > objective:
                    break
            return last

    cutoff = None
    if args.cutoff_eigen_ratio is not None:
        cutoff = CutoffEigenRatio(args.cutoff_eigen_ratio)
    elif args.cutoff_trace_ratio is not None:
        cutoff = CutoffTraceRatio(args.cutoff_trace_ratio)
    elif args.cutoff_threshold is not None:
        cutoff = CutoffThreshold(args.cutoff_threshold)
    else:
        raise Exceptions.InvalidArguments("Specify either cutoff_ratio or cutoff_threshold")
    return cutoff

def load_variance(path, trim_ensemble_version=True):
    logging.log(9, "Loading variance",)
    gene_variance_data = pandas.read_table(path)
    logging.log(9, "Indexing variance")
    if trim_ensemble_version:
        k = gene_variance_data.gene.str.split(".").str.get(0)
        gene_variance_data.gene = k
    gene_variance_data = gene_variance_data.set_index(["gene", "model"])
    gene_variance_data = gene_variance_data.sort_index()
    return gene_variance_data

def context_from_args(args):
    context = None

    logging.info("Loading genes")
    genes = PredictionModel.load_genes(args.models_folder)

    if args.model_product:
        logging.info("Context for model product")
        definition={
            MatrixManager.K_MODEL:"gene",
            MatrixManager.K_ID1:"model1",
            MatrixManager.K_ID2:"model2",
            MatrixManager.K_VALUE:"value"
        }
        matrix_manager = MatrixManager.load_matrix_manager(args.model_product, definition=definition, permissive=args.permissive_model_product)
        metaxcan_manager = MetaXcanResultsManager.build_manager(args.metaxcan_folder, filters=args.metaxcan_filter)
        cutoff = _cutoff(args)
        context = SimpleContext(metaxcan_manager, matrix_manager, genes, cutoff, args.regularization)
    elif args.snp_covariance:
        logging.info("Context for snp covariance")
        if args.cleared_snps:
            intersection = KeyedDataSource.load_data_column(args.cleared_snps, "rsid")
            intersection = set(intersection)
        else:
            intersection = GWASAndModels.gwas_model_intersection(args)

        if len(intersection) == 0:
            raise Exceptions.ReportableException("No intersection of snps between GWAS and models.")

        trim = not args.full_ensemble_id
        model_manager = PredictionModel.load_model_manager(args.models_folder, trim_ensemble_version=trim, Klass=PredictionModel._ModelManager)

        def _check_in(comps, intersection):
            return comps[1] not in intersection or comps[2] not in intersection

        snp_covariance_streamer = DataFrameStreamer.data_frame_streamer(args.snp_covariance, "GENE", additional_skip_row_check= lambda x: _check_in(x, intersection))

        cutoff = _cutoff(args)
        metaxcan_manager = MetaXcanResultsManager.build_manager(args.metaxcan_folder, filters=args.metaxcan_filter, file_name_pattern=args.metaxcan_file_name_parse_pattern)
        context = ExpressionStreamedContext(metaxcan_manager, snp_covariance_streamer, model_manager, genes, cutoff, args.regularization, args.full_ensemble_id)
        context.check()
    return context


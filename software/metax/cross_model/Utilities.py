import numpy
import pandas

from JointAnalysis import Context, ContextMixin
from .. import MatrixManager
from .. import Exceptions
from ..misc import DataFrameStreamer
from ..genotype import GenotypeAnalysis
from ..metaxcan import MetaXcanResultsManager

class SimpleContext(Context, ContextMixin):
    def __init__(self, metaxcan_results_manager, matrix_manager, cutoff, epsilon):
        self.metaxcan_results_manager = metaxcan_results_manager
        self.matrix_manager = matrix_manager
        self.model_name_index = _index(matrix_manager)
        # Yeah, a callable object. Go figure.
        self.cutoff = cutoff
        self.epsilon = epsilon

    def get_genes(self):
        matrix_genes = self.matrix_manager.model_labels()
        results_genes = self.metaxcan_results_manager.get_genes()
        genes = {x for x in matrix_genes if x.split(".")[0] in results_genes}
        return genes

class ExpressionStreamedContext(Context, ContextMixin):
    def __init__(self, metaxcan_results_manager, gene_variance_data, snp_covariance_streamer, cutoff, epsilon):
        self.metaxcan_results_manager = metaxcan_results_manager
        self.gene_variance_data = gene_variance_data,
        self.snp_covariance_streamer = snp_covariance_streamer
        self.cutoff = cutoff
        self.epsilon = epsilon
        self.matrix_manager = None

    def get_genes(self):
        for d in self.snp_covariance_streamer:
            self.matrix_manager = GenotypeAnalysis.GeneExpressionMatrixManager(self.gene_variance_data, d)
            yield d.GENE.values[0]

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

def context_from_args(args):
    context = None
    if args.model_product:
        definition={
            MatrixManager.K_MODEL:"gene",
            MatrixManager.K_ID1:"model1",
            MatrixManager.K_ID2:"model2",
            MatrixManager.K_VALUE:"value"
        }
        matrix_manager = MatrixManager.load_matrix_manager(args.model_product, definition=definition)
        metaxcan_manager = MetaXcanResultsManager.build_manager(args.metaxcan_folder, filters=args.metaxcan_filter)
        cutoff = _cutoff(args)
        context = SimpleContext(metaxcan_manager, matrix_manager, cutoff, args.regularization)
    elif args.expression_data_prefix:
        snp_covariance_path = args.expression_data_prefix + "_snp_covariance.txt.gz"
        snp_covariance_streamer = DataFrameStreamer.data_frame_streamer(snp_covariance_path, "GENE")

        gene_variance_path = args.expression_data_prefix + "_gene_variance.txt.gz"
        gene_variance_data = pandas.read_table(gene_variance_path)

        cutoff = _cutoff(args)
        metaxcan_manager = MetaXcanResultsManager.build_manager(args.metaxcan_folder, filters=args.metaxcan_filter)
        context = ExpressionStreamedContext(metaxcan_manager, gene_variance_data, snp_covariance_streamer, cutoff, args.regularization)
    return context


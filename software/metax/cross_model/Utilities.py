import numpy
from JointAnalysis import Context
from ..metaxcan import MetaXcanResultsManager
from .. import MatrixManager
from .. import Exceptions

class SimpleContext(Context):
    def __init__(self, metaxcan_results_manager, matrix_manager, cutoff):
        self.metaxcan_results_manager = metaxcan_results_manager
        self.matrix_manager = matrix_manager
        self.model_name_index = _index(matrix_manager)
        # Yeah, a callable object. Go figure.
        self.cutoff = cutoff

    def get_genes(self):
        matrix_genes = self.matrix_manager.model_labels()
        results_genes = self.metaxcan_results_manager.get_genes()
        genes = {x for x in matrix_genes if x.split(".")[0] in results_genes}
        return genes

    def get_metaxcan_zscores(self, gene):
        if "." in gene: gene = gene.split(".")[0]
        results = self.metaxcan_results_manager.results_for_gene(gene)
        return results

    def get_model_matrix(self, gene, tissues):
        return self.matrix_manager.get(gene, tissues)

    def get_cutoff(self, matrix):
        return self.cutoff(matrix)

def _index(matrix_manager):
    labels = matrix_manager.model_labels()
    labels = {x.split(".")[0]:x for x in labels}
    return labels

def _cutoff(args):
    class CutoffRatio(object):
        def __init__(self, cutoff_ratio):
            self.cutoff_ratio = float(cutoff_ratio)

        def __call__(self, matrix):
            trace = numpy.trace(matrix)
            cutoff = self.cutoff_ratio * trace
            return cutoff

    class CutoffThreshold(object):
        def __init__(self, cutoff_threshold):
            self.cutoff_threshold = float(cutoff_threshold)

        def __call__(self, matrix):
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
    if args.cutoff_ratio:
        cutoff = CutoffRatio(args.cutoff_ratio)
    elif args.cutoff_threshold:
        cutoff = CutoffThreshold(args.cutoff_threshold)
    else:
        raise Exceptions.InvalidArguments("Specify either cutoff_ratio or cutoff_threshold")
    return cutoff

def context_from_args(args):
    definition={
        MatrixManager.K_MODEL:"gene",
        MatrixManager.K_ID1:"model1",
        MatrixManager.K_ID2:"model2",
        MatrixManager.K_VALUE:"value"
    }
    matrix_manager = MatrixManager.load_matrix_manager(args.model_product, definition=definition)
    metaxcan_manager = MetaXcanResultsManager.build_manager(args.metaxcan_folder, filters=args.metaxcan_filter)
    cutoff = _cutoff(args)
    context = SimpleContext(metaxcan_manager, matrix_manager, cutoff)
    return context


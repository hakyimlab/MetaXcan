from JointAnalysis import Context
from ..metaxcan import MetaXcanResultsManager
from .. import MatrixManager

class SimpleContext(Context):
    def __init__(self, metaxcan_results_manager, matrix_manager):
        self.metaxcan_results_manager = metaxcan_results_manager
        self.matrix_manager = matrix_manager
        self.model_name_index = _index(matrix_manager)

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

def _index(matrix_manager):
    labels = matrix_manager.model_labels()
    labels = {x.split(".")[0]:x for x in labels}
    return labels

def context_from_args(args):
    definition={
        MatrixManager.K_MODEL:"gene",
        MatrixManager.K_ID1:"model1",
        MatrixManager.K_ID2:"model2",
        MatrixManager.K_VALUE:"value"
    }
    matrix_manager = MatrixManager.load_matrix_manager(args.model_product, definition=definition)
    metaxcan_manager = MetaXcanResultsManager.build_manager(args.metaxcan_folder, filters=args.metaxcan_filter)
    context = SimpleContext(metaxcan_manager, matrix_manager)
    return context


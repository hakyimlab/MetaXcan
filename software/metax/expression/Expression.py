from .. import Exceptions

class ExpressionManager(object):
    def expression_for_gene(self, gene): raise Exceptions.ReportableException("Unimplemented expression_for_gene")
    def get_genes(self): raise Exceptions.ReportableException("Unimplemented get_genes")
    def enter(self): pass
    def exit(self): pass
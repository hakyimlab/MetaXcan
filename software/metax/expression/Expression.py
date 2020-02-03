from .. import Exceptions

class ExpressionManager(object):
    def expression_for_gene(self, gene): raise Exceptions.ReportableException("Unimplemented expression_for_gene")
    def get_genes(self): raise Exceptions.ReportableException("Unimplemented get_genes")
    def enter(self): pass
    def exit(self): pass

class Expression(object):
    def expression_for_gene(self, gene): raise Exceptions.ReportableException("Unimplemented expression_for_gene")
    def get_genes(self): raise Exceptions.ReportableException("Unimplemented get_genes")
    def enter(self): pass
    def exit(self): pass

class ExpressionFromData(Expression):
    def __init__(self, data):
        self.d = data

    def expression_for_gene(self, gene):
        k = self.d[gene]
        return k

    def get_genes(self):
        return list(self.d.keys())

    def enter(self):
        pass

    def exit(self):
        pass

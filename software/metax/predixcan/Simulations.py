__author__ = "alvaro barbeira"

import numpy
import pandas

from . import Utilities
from . import MultiPrediXcanAssociation
from ..misc import Math
from ..expression import HDF5Expression

def mp_callback(model, result, vt_projection, model_keys, save):
    coefs = MultiPrediXcanAssociation._coefs(result, vt_projection, model_keys)
    save["coefs"] = coefs

class Context(object):
    def __init__(self, expression_manager, phenotype_generator, filter):
        self.expression_manager = expression_manager
        self.phenotype_generator = phenotype_generator
        self.filter = filter

    def get_genes(self):
        return self.expression_manager.get_genes()[0:10]

    def __enter__(self):
        self.expression_manager.enter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.expression_manager.exit()

    def get_mp_simulation(self, gene):
        expression = self.expression_manager.expression_for_gene(gene)
        phenotype, description = self.phenotype_generator.get(expression, gene)
        if phenotype is None:
            return None, None
        return Utilities.DumbMTPContext(expression, phenotype, gene, self.filter), description

class PhenotypeGenerator(object):
    def __init__(self): raise RuntimeError("Not implemented")

class RandomPhenotypeGenerator(PhenotypeGenerator):
    def __init__(self):
        pass

    def get(self, expression, gene):
        k = list(expression.keys())[0]
        e = expression[k]
        n  = len(e)
        pheno = numpy.random.uniform(size=n)
        description = pandas.DataFrame({ "variable":["covariate"], "param": [1.0]})
        return pheno, description


class LinearCombinationPhenotypeGenerator(PhenotypeGenerator):
    def __init__(self, combination):
        self.combination = combination

    def get(self, expression, gene):
        return _pheno_from_combination(expression, self.combination)

def _pheno_from_combination(expression, combination):
    ok = True
    _k = list(expression.keys())[0]
    _e = expression[_k]
    n = len(_e)
    e = numpy.zeros(n)
    used = set()
    for k, v in combination.iteritems():
        if k in expression:
            e +=  expression[k] * v
            used.add(k)
        elif "covariate" in k:
            e += numpy.random.normal(size=n) * v
            used.add(k)
        else:
            # If we couldn't build a model with the desired combination, abort
            ok = False
            break

    if not ok:
        return None, None

    _c = {x: v for x, v in combination.iteritems() if x in used}
    pheno = e
    description = pandas.DataFrame({"variable": list(_c.keys()), "param": list(_c.values())})
    return pheno, description

def context_from_args(args):
    expression = HDF5Expression.ExpressionManager(args.hdf5_expression_folder, args.expression_pattern, code_999=args.code_999, standardise= args.standardize_expression)
    if args.simulation_type == "random":
        phenotype = RandomPhenotypeGenerator()
    elif args.simulation_type == "combination":
        _c = {"Adipose_Subcutaneous":1.0, "Brain_Cerebellum":1.0, "covariate":1.0}
        phenotype = LinearCombinationPhenotypeGenerator(_c)
    else:
        raise RuntimeError("Wrong phenotype simulation spec")
    filter = Utilities._filter_from_args(args)
    context = Context(expression, phenotype, filter)
    return context


def simulate(gene, context):
    save_results = {}
    _cb = lambda model, result, vt_projection, model_keys: mp_callback(model, result, vt_projection, model_keys, save_results)
    _context, _description = context.get_mp_simulation(gene)
    if _context is None:
        return None, None

    r = MultiPrediXcanAssociation.multi_predixcan_association(gene, _context, _cb)
    description = _description.assign(gene=gene, type="truth")
    coefs = save_results["coefs"].assign(gene=gene, type="result")
    description = pandas.concat([description, coefs])
    return r, description
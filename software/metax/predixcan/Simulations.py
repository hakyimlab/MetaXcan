__author__ = "alvaro barbeira"

import numpy
import pandas
import copy
import math

from . import Utilities
from . import MultiPrediXcanAssociation, PrediXcanAssociation
from ..expression import HDF5Expression, Expression

########################################################################################################################
def mp_callback(gene, model, result, vt_projection, variance, model_keys, coefs, save):
    save["coefs"] = coefs

class Context(object):
    def __init__(self, expression_manager, phenotype_generator, filter, do_predixcan=False, only_truth=False):
        self.expression_manager = expression_manager
        self.phenotype_generator = phenotype_generator
        self.filter = filter
        self.do_predixcan = do_predixcan
        self.only_truth = only_truth

    def do_predixcan(self):
        return self.do_predixcan

    def get_genes(self):
        return self.expression_manager.get_genes()

    def __enter__(self):
        self.expression_manager.enter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.expression_manager.exit()

    def get_mp_simulation(self, gene):
        if not gene:
            return Utilities.DumbMTPContext(None, None, None, self.filter), None, None

        expression = self.expression_manager.expression_for_gene(gene)
        phenotype, description = self.phenotype_generator.get(expression, gene)
        if phenotype is None:
            return None, None, None

        _cp = None
        if self.do_predixcan:
            _cp = {}
            for t in description.itertuples():
                if "covariate" in t.variable: continue
                _cp[t.variable] = Utilities.DumbPContext(expression[t.variable], phenotype, gene, self.filter)

        if self.only_truth:
            expression = {x.variable:expression[x.variable] for x in description.itertuples() if x.variable in expression}

        return Utilities.DumbMTPContext(expression, phenotype, gene, self.filter), _cp, description

########################################################################################################################

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
    def __init__(self, combination, covariate_sd, use_all=None):
        self.combination = combination
        self.covariate_sd = covariate_sd
        self.use_all = use_all

    def get(self, expression, gene):
        combination = copy.copy(self.combination)
        if self.use_all:
            if type(self.use_all) == float:
                c = self.use_all
            elif self.use_all == "ONE_VAR":
                c = math.sqrt(1.0 / len(expression))
            elif self.use_all == "FIX_VAR":
                c = 1.0
                combination["covariate"] = math.sqrt(len(expression)*99)
            else:
                raise RuntimeError("Unsupported option")
            for e in list(expression.keys()):
                combination[e] = c

        return _pheno_from_combination(expression, combination, self.covariate_sd)

class CombinationOfCorrelatedPhenotypeGenerator(PhenotypeGenerator):
    def __init__(self, covariate_coefficient, covariate_sd, threshold):
        self.covariate_coefficient = covariate_coefficient
        self.threshold = threshold
        self.covariate_sd = covariate_sd

    def get(self, expression, gene):
        # Get the tissue with the most correlated siblings;
        # then average them build a phenotype
        values = list(expression.values())
        if len(values) == 1:
            return None, None
        e = values
        c = numpy.corrcoef(e)
        d = len(expression)
        f = 0
        r = 0
        for i in range(0, d):
            f_ = numpy.sum(c[i] > self.threshold)
            if f_ > f:
                r = i
                f = f_

        if f<2:
            return None, None

        which = c[r] > self.threshold
        keys = list(expression.keys())
        combination = {keys[i]:math.sqrt(1.0/f) for i in range(0, d) if which[i]}
        #for i in xrange(0,d):
        #    combination["covariate_{}".format(i)] = 10.0/f
        combination["covariate"] = self.covariate_coefficient

        return  _pheno_from_combination(expression, combination, self.covariate_sd)

def _pheno_from_combination(expression, combination, covariate_sd):
    ok = True
    _k = list(expression.keys())[0]
    _e = expression[_k]
    n = len(_e)
    e = numpy.zeros(n)
    used = set()
    for k, v in combination.items():
        if k in expression:
            e +=  expression[k] * v
            used.add(k)
        elif "covariate" in k:
            e += numpy.random.normal(scale=covariate_sd, size=n) * v
            used.add(k)
        else:
            # If we couldn't build a model with the desired combination, abort
            ok = False
            break

    if not ok:
        return None, None

    _c = {x: v for x, v in combination.items() if x in used}
    pheno = e
    description = pandas.DataFrame({"variable": list(_c.keys()), "param": list(_c.values())})
    return pheno, description

########################################################################################################################

class SExpressionManager(Expression.ExpressionManager):
    def __init__(self, em):
        self.em = em
        self.which = None

    def expression_for_gene(self, gene):
        e = self.em.expression_for_gene(gene)

        if self.which is None:
            n = len(e[list(e.keys())[0]])
            s = 10000
            #self.which = numpy.random.choice([True, False], size=n, p=[s*1.0/n, 1 - s*1.0/n])
            self.which = list(numpy.random.choice(range(0,n), size=s, replace=False))

        e = {k:v[self.which] for k,v in e.items()}
        return e

    def get_genes(self): return self.em.get_genes()
    def enter(self): return self.em.enter()
    def exit(self): self.em.exit()


def context_from_args(args):
    #expression_ = HDF5Expression.ExpressionManager(args.hdf5_expression_folder, args.expression_pattern, code_999=args.code_999, standardise= args.standardize_expression)
    #expression = SExpressionManager(expression_)
    expression = HDF5Expression.ExpressionManager(args.expression_folder, args.expression_pattern,
                                                   code_999=args.code_999, standardise=args.standardize_expression)

    def _argumentize(x, t, default=1.0):
        return t(x) if x is not None else default

    p = {x[0]: x[1] for x in args.simulation_parameters} if args.simulation_parameters else {}
    covariate_coefficient = _argumentize(p.get("covariate_coefficient"), float)
    covariate_sd = _argumentize(p.get("covariate_sd"), float)
    if args.simulation_type == "random":
        phenotype = RandomPhenotypeGenerator()
    elif args.simulation_type == "combination":
        use_all = None
        if "model_spec" in p:
            _c = p.get("model_spec")
            if not _c:
                _c = {}
            else:
                _c = _c.split()
                _c = {_c[i*2]:float(_c[i*2+1]) for i in range(0,len(_c)/2)}
        elif "use_tissues" in p:
            _c = p.get("use_tissues").strip().split()
            _c = {x:math.sqrt(1.0/len(_c)) for x in _c}
        elif "use_all" in p:
            _c = {}
            if p["use_all"] == "ONE_VAR" or p["use_all"] == "FIX_VAR":
                use_all = p["use_all"]
            else:
                use_all = float(p["use_all"])
        _c["covariate"] = covariate_coefficient
        phenotype = LinearCombinationPhenotypeGenerator(_c, covariate_sd=covariate_sd, use_all=use_all)
    elif args.simulation_type == "combination_from_correlated":
        threshold = _argumentize(p.get("threshold"), float, 0.9)
        phenotype = CombinationOfCorrelatedPhenotypeGenerator(covariate_coefficient=covariate_coefficient, covariate_sd=covariate_sd, threshold=threshold)
    else:
        raise RuntimeError("Wrong phenotype simulation spec")
    filter = Utilities._filter_from_args(args)
    context = Context(expression, phenotype, filter, args.do_predixcan, args.only_truth)
    return context


########################################################################################################################

def simulate(gene, context):
    save_results = {}
    _cb = lambda gene, model, result, vt_projection, variance, model_keys, coefs: mp_callback(gene, model, result, vt_projection, variance, model_keys, coefs, save_results)

    _context_mt, _context_p, _description = context.get_mp_simulation(gene)
    if _context_mt is None:
        return None, None, None

    p = None
    if _context_p:
        p = pandas.DataFrame()
        for model,_c in _context_p.items():
            p_ = PrediXcanAssociation.predixcan_association(gene, _c)
            p_ = PrediXcanAssociation.dataframe_from_results([p_])
            p_["model"] = model
            p = pandas.concat([p,p_])


    r = MultiPrediXcanAssociation.multi_predixcan_association(gene, _context_mt, [_cb])
    description = _description.assign(gene=gene, type="truth")
    coefs = save_results["coefs"].assign(gene=gene, type="result")
    description = pandas.concat([description, coefs])

    return r, description, p
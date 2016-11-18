__author__ = 'heroico'

import logging
import numpy
import math
import Exceptions
import KeyedDataSet

BETA_Z = "beta_z"
BETA_Z_SIGMA_REF = "beta_z_and_ref"
METAXCAN = "metaxcan"
METAXCAN_FROM_REFERENCE = "metaxcan_from_reference"

def ZScoreScheme(name):
    scheme = None
    if name == BETA_Z:
        scheme = _BetaZ()
    elif name == BETA_Z_SIGMA_REF:
        scheme = _BetaZAndRef()
    elif name == METAXCAN:
        scheme = _MetaXcan()
    elif name == METAXCAN_FROM_REFERENCE:
        scheme = _MetaXcanFromReference()
    else:
        raise Exception("Unknown zscore scheme %s", name if name else "None")
    return scheme

class ZScoreCalculation(object):
    def __call__(self, gene, weights, beta_sets, covariance_matrix, valid_rsids):
        raise Exceptions.ReportableException("Wrong Zscore calculation scheme!")
        return  None, None, None

    def getValue(self, set, rsid):
        if not rsid in set.values_by_key:
            logging.log(5, "rsid %s not in %s", rsid, "set" if set.name is None else set.name)
            return None
        value = set.values_by_key[rsid]
        if value == "NA":
            logging.log(5, "rsid %s doesnt have a value for %s", rsid, "set" if set.name is None else set.name)
            return None
        value = float(value)
        return value

    def get_beta_z(self, beta_sets, rsid):
        b_z = beta_sets["beta_z"] if "beta_z" in beta_sets else KeyedDataSet.EMPTY
        value = self.getValue(b_z, rsid)
        return value

    def get_beta(self, beta_sets, rsid):
        betas = beta_sets["beta"] if "beta" in beta_sets else KeyedDataSet.EMPTY
        value = self.getValue(betas, rsid)
        return value

    def get_sigma_l(self, beta_sets, rsid):
        sigma_l = beta_sets["sigma_l"] if "sigma_l" in beta_sets else KeyedDataSet.EMPTY
        value = self.getValue(sigma_l, rsid)
        return value

    def get_reference_sigma_l(self, variances, rsid):
        if not rsid in variances:
            logging.log(5, "rsid %s not in variances", rsid)
            return None

        sigma_l = math.sqrt(variances[rsid])
        return sigma_l

class _BetaZ(ZScoreCalculation):
    def __call__(self, gene, weights, beta_sets, covariance_matrix, valid_rsids):
        #TODO: have beta_z_validation be provided by each calculation class,
        #So as to check for standard error in input sets if applicable
        weight_values, variances = preProcess(covariance_matrix, valid_rsids, weights, beta_sets, beta_z_validation)

        effect_size="NA"
        zscore = "NA"
        n = 0
        #dot_product is Var(g)
        dot_product = numpy.dot(numpy.dot(numpy.transpose(weight_values), covariance_matrix), weight_values)
        if dot_product > 0:
            denominator = math.sqrt(float(dot_product))
            zscore_numerator_terms = []
            effect_size_numerator_terms = []
            for rsid in valid_rsids:
                w = weights[rsid].weight

                b_z = self.beta_z(beta_sets, rsid)
                if b_z is None:
                    continue

                s_l = self.sigma_l(beta_sets, variances, rsid)
                if not s_l:
                    continue

                z_term = w * b_z * s_l
                zscore_numerator_terms.append(z_term)

                #effect size
                b = self.get_beta(beta_sets, rsid)
                if not b:
                    continue
                e_term = w * b * s_l**2
                effect_size_numerator_terms.append(e_term)

            n = len(zscore_numerator_terms)
            if n > 0:
                numerator = sum(zscore_numerator_terms)
                zscore = str(numerator/denominator)
                if len(effect_size_numerator_terms):
                    e_numerator = sum(effect_size_numerator_terms)
                    effect_size = str(e_numerator/dot_product)
            else:
                logging.log(7,"No terms for %s ", gene)

        return zscore, str(n), str(dot_product), effect_size

    def beta_z(self, beta_sets, rsid):
        return self.get_beta_z(beta_sets, rsid)

    def sigma_l(self, beta_sets, variances, rsid):
        return self.get_sigma_l(beta_sets, rsid)

class _BetaZAndRef(_BetaZ):
    def sigma_l(self, beta_sets, variances, rsid):
        return self.get_reference_sigma_l(variances, rsid)

class _MetaXcan(ZScoreCalculation):
    def __call__(self, gene, weights, beta_sets, covariance_matrix, valid_rsids):
        #TODO: have beta_z_validation be provided by each calculation class,
        #So as to check for standard error in input sets if applicable
        weight_values, variances = preProcess(covariance_matrix, valid_rsids, weights, beta_sets, beta_validation)

        pre_zscore = "NA"
        effect_size = "NA" # Stubbed out, proper approach in this scheme is not clear
        n = 0
        #dot_product is Var(g)
        dot_product = numpy.dot(numpy.dot(numpy.transpose(weight_values), covariance_matrix), weight_values)
        if dot_product > 0:
            denominator = math.sqrt(float(dot_product))
            numerator_terms = []
            for rsid in valid_rsids:
                w = weights[rsid].weight

                b = self.get_beta(beta_sets, rsid)
                if b is None:
                    continue

                s_l = self.sigma_l(beta_sets, variances, rsid)
                if not s_l:
                    continue

                term = w * b * s_l**2
                numerator_terms.append(term)

            n = len(numerator_terms)
            if n > 0:
                numerator = sum(numerator_terms)
                pre_zscore = str(numerator/denominator)
            else:
                logging.log(7, "No terms for %s ", gene)

        return pre_zscore, str(n), str(dot_product), effect_size

    def sigma_l(self, beta_sets, variances, rsid):
        return self.get_sigma_l(beta_sets, rsid)

class _MetaXcanFromReference(_MetaXcan):
    def sigma_l(self, beta_sets, variances, rsid):
        return self.get_reference_sigma_l(variances, rsid)

def check_input_set_rsid(input_set, rsid):
    if not rsid in input_set.values_by_key or input_set.values_by_key[rsid] == "NA":
        return False
    return True

def beta_validation(beta_sets, rsid):
    if not "beta" in beta_sets:
        return False

    b = beta_sets["beta"]
    if not check_input_set_rsid(b, rsid):
        return False

    return True

def beta_z_validation(beta_sets, rsid):
    if not "beta_z" in beta_sets:
        return False

    b_z = beta_sets["beta_z"]
    if not check_input_set_rsid(b_z, rsid):
        return False

    return True

def preProcess(covariance_matrix, valid_rsids, weights, beta_sets, validation):
    weight_values = []
    variances = {}
    for i,rsid in enumerate(valid_rsids):
        if rsid not in weights:
            raise Exceptions.ReportableException("RSID %s can't be found in the weights database. Are you sure your covariance data matches the weights database you are using?" % (rsid))
        weight = weights[rsid].weight
        if not validation(beta_sets, rsid):
            logging.log(7, "snp %s not present in beta data, skipping weight")
            # this will effectively skip this rsid at (w * G * w)
            weight = 0
        weight_values.append(weight)
        if covariance_matrix.ndim == 0:
            variances[rsid] = float(covariance_matrix)
        else:
            variances[rsid] = covariance_matrix[i][i]
    return weight_values, variances

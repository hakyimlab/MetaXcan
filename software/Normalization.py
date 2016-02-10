__author__ = 'heroico'

import numpy
import logging
import math

NONE="none"
FROM_PHENO="from_pheno"
FROM_REFERENCE="from_reference"

def normalizationScheme(scheme, covariances=None, weight_db_logic=None):
    n = None
    if scheme == NONE:
        n = NoNormalization()
    elif scheme == FROM_PHENO:
        n = _BetaNormalization()
    elif scheme == FROM_REFERENCE:
        n = _ReferenceNormalization(covariances, weight_db_logic)
    else:
        raise Exception("Unknown normalization: %s", scheme if scheme else None)
    return n

class Normalization(object):
    def __init__(self):
        pass

    def update(self, beta_sets):
        pass

    def calculateNormalization(self):
        return 1

class NoNormalization(Normalization):
    pass

class _BetaNormalization(Normalization):
    def __init__(self):
        self.ses = []
        self.sigmas = []

    def update(self, beta_sets):
        se = beta_sets["se"]
        self.ses.append(se)

        sigma = beta_sets["sigma_l"]
        self.sigmas.append(sigma)

    def calculateNormalization(self):
        logging.info("Calculating normalization from phenotype")
        y = []
        x = []
        for i, ses_data in enumerate(self.ses):
            logging.log(6,"processing standard error %i", i)
            sigma_data = self.sigmas[i]
            for j in xrange(0, len(ses_data.data)):
                se = ses_data.data[j]
                sigma = sigma_data.data[j]
                if se == "NA" or sigma == "NA":
                    continue

                y.append(1/float(se))
                s = float(sigma)
                x.append(s)

        x = numpy.array(x)
        x = x[:,numpy.newaxis]
        y = numpy.array(y)
        a = numpy.linalg.lstsq(x,y)[0]
        return float(a)

class _ReferenceNormalization(Normalization):
    def __init__(self, covariances, weight_db_logic):
        self.ses = []
        self.covariances = covariances
        self.weight_db_logic = weight_db_logic

    def update(self, beta_sets):
        se = beta_sets["se"]
        self.ses.append(se)

    def calculateNormalization(self):
        logging.info("Calculating normalization from reference")
        y = []
        x = []
        for i, ses_data in enumerate(self.ses):
            for j in xrange(0, len(ses_data.data)):
                se = ses_data.data[j]
                if se == "NA":
                    continue

                rsid = ses_data.keys[j]
                genes = self.weight_db_logic.genes_for_an_rsid[rsid]
                gene_count = len(genes)
                if not gene_count:
                    logging.log(5, "no genes for rsid %s, skipping", rsid)
                    continue

                entry = None
                for gene in genes:
                    if gene in self.covariances:
                        entry = self.covariances[gene]
                        logging.log(5, "picked gene %s for rsid %s from %d", gene, rsid, gene_count)
                        break

                if not entry:
                    logging.log(8, "rsid %s without proper covariance matrix, skipping", rsid)
                    continue

                covariance_matrix = entry[0]
                valid_rsids = entry[1]
                if not rsid in valid_rsids:
                    logging.log(8, "rsid %s not in covariance matrix, skipping", rsid)
                    continue

                index = valid_rsids.index(rsid)
                sigma = math.sqrt(covariance_matrix[index][index])
                s = float(sigma)
                x.append(s)

                y.append(1/float(se))

        x = numpy.array(x)
        x = x[:,numpy.newaxis]
        y = numpy.array(y)
        a = numpy.linalg.lstsq(x,y)[0]
        return float(a)

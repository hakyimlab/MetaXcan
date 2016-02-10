#! /usr/bin/env python
__author__ = 'heroico'

import os
import logging
import Logging
import Utilities
import MatrixUtilities
import KeyedDataSet
import WeightDBUtilities
import GWASUtilities

if __name__ == "__main__":
    Logging.configureLogging(logging.INFO)

COV = "intermediate/cov/covariance.txt.gz"
BETA = "data/T1D-GWAS"
WEIGHT_DB_PATH = "data/DGN-WB_0.5.db"

def loadDosageFile(path):
    callback = GWASUtilities.GWASSNPInfoLineCollector()
    dosage_loader = GWASUtilities.GWASDosageFileLoader(path, True, callback)
    keyed_data_set = dosage_loader.load()
    return keyed_data_set

logging.info("Loading weight db")
weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(WEIGHT_DB_PATH)

logging.info("Loading covariance file")
covariance_contents = MatrixUtilities.loadMatrixFromFile(COV)

logging.info("Loading betas")
beta_contents = Utilities.contentsWithPatternsFromFolder(BETA, [".gz"])
results = []
for beta_name in beta_contents:
    logging.info("Processing %s", beta_name)
    beta_path = os.path.join(BETA, beta_name)

    beta_data = loadDosageFile(beta_path)[0]
    for snp, value in beta_data.values_by_key.iteritems():
        if not snp in weight_db_logic.genes_for_an_rsid:
            logging.log(7, "rsid %s not found in DB", snp)
            continue
        genes = weight_db_logic.genes_for_an_rsid[snp]
        if not genes:
            logging.info("no gene for %s", snp)
            continue
        gene = genes[0]
        if not gene in covariance_contents:
            logging.info("(%s, %s) without covariance matrix", gene, snp)
            continue

        covariance_matrix, valid_rsids = covariance_contents[gene]
        if not snp in valid_rsids:
            logging.log(9, "rsid %s not in matrix", snp)
            continue
        i = valid_rsids.index(snp)
        var = float(covariance_matrix[i][i])

        freq = float(value[GWASUtilities.GWASSNPInfoLineCollector.FREQ])
        s2 = 2 * freq * (1 - freq)
        results.append((snp, var, s2))

with open("results/kk.csv", "wb") as file:
    file.write("rsid,ref_var,freq_var\n")
    for result in results:
        snp = result[0]
        var = result[1]
        s2 = result[2]
        line = "%s,%f,%f\n" % (snp, var, s2)
        file.write(line)

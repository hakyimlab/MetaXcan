#! /usr/bin/env python
__author__ = 'heroico'
import logging
import os
import numpy
import WeightDBUtilities
import Logging
import Utilities
import MatrixUtilities
import ZScoreCalculation


class GeneStats(object):
    def __init__(self, args):
        self.weight_db_path = args.weight_db_path
        self.folder_covariance = args.covariance_folder
        self.output_file = args.output_file

    def run(self):
        logging.info("Loading weight db")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.weight_db_path)

        logging.info("Loading covariance file")
        file = Utilities.contentsWithPatternsFromFolder(self.folder_covariance, [".gz"])[0]
        path = os.path.join(self.folder_covariance, file)
        covariance_contents = MatrixUtilities.loadMatrixFromFile(path)

        logging.info("Getting stats")
        results = []
        for gene, entry in covariance_contents.iteritems():
            covariance_matrix = entry[0]
            valid_rsids = entry[1]

            weights = weight_db_logic.weights_by_gene_name[gene]
            weight_values, variances = ZScoreCalculation.preProcess(covariance_matrix, valid_rsids, weights)

            w_w = numpy.dot(numpy.transpose(weight_values), weight_values)
            dot_product = numpy.dot(numpy.dot(numpy.transpose(weight_values), covariance_matrix), weight_values)
            det = numpy.linalg.det(covariance_matrix)

            eigenvalues, eigenvectors = numpy.linalg.eigh(covariance_matrix)
            eigenmax = numpy.amax(eigenvalues)
            eigenmin = numpy.amin(eigenvalues)
            n_small = 0
            for eigen in eigenvalues:
                if eigen < 1e-7:
                    n_small += 1
            diag = covariance_matrix.diagonal()
            mean_var = numpy.mean(diag)

            line = (gene, str(len(weight_values)), str(float(dot_product)), str(float(det)), str(float(w_w)), str(float(mean_var)), str(float(eigenmin)), str(float(eigenmax)), str(n_small))
            results.append(line)

#gene, n.snps, WW, W\Gamma W, eig(\Gamma).max, eig(\Gamma).min, #eigs<1e-8, VAR_g, zscore_g
        logging.info("saving results")
        with open(self.output_file, "w") as file:
            header = ",".join(["gene", "m_snp_count", "w_gamma_w", "det", "w_w", "mean_var", "eigenmin", "eigenmax", "n_eigen_e-7"])+"\n"
            file.write(header)
            for line in results:
                text = ",".join(line)+"\n"
                file.write(text)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build betas from GWAS data.')

    parser.add_argument("--weight_db_path",
                        help="name of weight db in data folder",
                        default="data/cross-tissue_0.5.db")

    parser.add_argument("--covariance_folder",
                        help="name of folder containing covariance data",
                        default="intermediate/COV_CROSSTISSUE")

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="forensics_ct/gene_stats.csv")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = GeneStats(args)
    work.run()


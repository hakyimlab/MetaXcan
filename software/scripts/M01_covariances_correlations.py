#!/usr/bin/env python
__author__ = 'heroico'

import gc
import logging
import numpy
import sqlite3
import os
import gzip
import ntpath
import metax.WeightDBUtilities as WeightDBUtilities
import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.ThousandGenomesUtilities as ThousandGenomesUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.Formats as Formats


def pathLeaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

class ProcessWeightDB(object):
    def __init__(self, args):
        self.weight_db = pathLeaf(args.weight_db)
        self.db_path = args.weight_db
        self.data_folder = args.input_folder
        self.correlation_output_folder = args.correlation_output_folder
        self.output_correlations = args.correlations
        self.covariance_output_folder = args.covariance_output_folder
        self.output_covariances = args.covariances

        self.input_format = args.input_format
        self.output_format = args.output_format

        self.found_genes_for_covariance = {}
        self.found_genes_for_correlation = {}

        self.min_maf_filter = float(args.min_maf_filter) if args.min_maf_filter else None
        self.max_maf_filter = float(args.max_maf_filter) if args.max_maf_filter else None

    def run(self):
        if not self.output_correlations and not self.output_covariances:
            return

        logging.info("Loading Weights")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.db_path)

        if not os.path.exists(self.correlation_output_folder) and self.output_correlations:
            os.makedirs(self.correlation_output_folder)

        if not os.path.exists(self.covariance_output_folder) and self.output_covariances:
            os.makedirs(self.covariance_output_folder)

        if self.output_format == Formats.DatabaseFile:
            self.buildDBS(weight_db_logic)
        elif self.output_format == Formats.FlatFile:
            self.buildFiles(weight_db_logic)
        else:
            logging.info("Wrong output format: %s", self.output_format)

    def buildDBS(self, weight_db_logic):
        names = Utilities.dosageNamesFromFolder(self.data_folder)
        for name in names:
            snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)

            if self.output_correlations:
                self.buildCorrelationDBS(weight_db_logic, name, snps, snps_by_rsid)

            if self.output_covariances:
                self.buildCovarianceDBS(weight_db_logic, name, snps, snps_by_rsid)

            gc.collect() #Python may have been saving numbers for reuse. Pretty please ask him to delete it. Files might be humongous.

    def buildFiles(self, weight_db_logic):
        do_correlations = self.output_correlations
        if self.output_correlations:
            correlation_path = self.correlationFilePath()
            if os.path.exists(correlation_path):
                logging.info("%s already exists, delete it if you want it figured out again", correlation_path)
                do_correlations = False
            else:
                self.writeFileHeader(correlation_path)

        do_covariances = self.output_covariances
        if self.output_covariances:
            covariance_path = self.covarianceFilePath()
            if os.path.exists(covariance_path):
                logging.info("%s already exists, delete it if you want it figured out again", covariance_path)
                do_covariances = False
            else:
                self.writeFileHeader(covariance_path)

        if not do_covariances and not do_correlations:
            return

        names = Utilities.dosageNamesFromFolder(self.data_folder)
        for name in names:
            snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)
            if do_correlations and self.output_correlations:
                self.addToCorrelationFile(weight_db_logic, name, snps, snps_by_rsid)

            if do_covariances and self.output_covariances:
                self.addToCovarianceFile(weight_db_logic, name, snps, snps_by_rsid)

    def writeFileHeader(self,path):
        with gzip.open(path, "ab") as file:
            file.write("GENE RSID1 RSID2 VALUE\n")

    def correlationFilePath(self):
        path = os.path.join(self.correlation_output_folder, "correlation.txt.gz")
        return path

    def covarianceFilePath(self):
        path = os.path.join(self.covariance_output_folder, "covariance.txt.gz")
        return path

    def getSNPS(self, name, weight_db_logic):
        dosageLoader = None
        if self.input_format == Formats.IMPUTE:
            dosageLoader = ThousandGenomesUtilities.IMPUTEDosageLoader(self.data_folder, name) #outdated code
        elif self.input_format == Formats.PrediXcan:
            dosageName = Utilities.dosageName(name)
            path = os.path.join(self.data_folder, dosageName)
            dosageLoader = PrediXcanFormatUtilities.PrediXcanFormatDosageLoader(path, weight_db_logic)
        else:
            logging.info("Invalid input format: %s", self.input_format)
            return
        snps, snps_by_rsid = dosageLoader.load()
        return snps, snps_by_rsid

    def addToCovarianceFile(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Adding to covariance for %s-%s", name, self.weight_db)

        genes = weight_db_logic.weights_by_gene_name.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            if percent == last_reported_percent+10:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent

            entries = self.buildCovarianceEntries(name, gene, weight_db_logic, snps_by_rsid)

            if len(entries) == 0:
                logging.log(6,"Gene %s has no snps in current file", gene)
                continue

            covariance_path = self.covarianceFilePath()
            self.addToFile(covariance_path, gene, entries)

    def addToFile(self, path, gene, entries):
        with gzip.open(path, "ab") as file:
            for entry in entries:
                line = " ".join([gene, entry[0], entry[1], entry[2]])+"\n"
                file.write(line)

    def buildCovarianceDBS(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Building covariance databases for %s-%s", name, self.weight_db)

        genes = weight_db_logic.weights_by_gene_name.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            if percent == last_reported_percent+10:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent

            output_name_base = "cov-"+name+"-"+gene+"-"+self.weight_db
            output_name = self.covariance_output_folder + "/" + output_name_base
            if os.path.exists(output_name):
                logging.info("Covariance matrix for %s already exists, delete if you want it figured out again", gene)
                continue

            entries = self.buildCovarianceEntries(name, gene, weight_db_logic, snps_by_rsid)

            if len(entries) == 0:
                logging.log(6,"Gene %s has no snps in current file", gene)
                continue

            self.writeCovarianceDB(output_name, entries)

    def buildCovarianceEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        weights_in_gene = weight_db_logic.weights_by_gene_name[gene]
        rsids_from_genes = weights_in_gene.keys()

        #gather as much data as we can work on
        related_rsids, related_data = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)

        if len(related_rsids) == 0:
            return []

        self.updateFoundCovariance(gene, name)

        #covariance matrix of related SNP's data
        array = numpy.array(related_data)
        cov = numpy.cov(array)

        #translate into sql entries
        entries = self.buildMatrixOutputEntries(cov, rsids_from_genes, related_rsids, snps_by_rsid)
        if not len(entries):
            raise NameError("Couldn not build covariance entries for (%s,%s)" %(name,gene))
        return entries

    def updateFoundCovariance(self, gene, name):
        found = None
        if gene in self.found_genes_for_covariance:
            found = self.found_genes_for_covariance[gene]
            logging.info("Gene %s found again for %s", gene, name)
        else:
            found = []
            self.found_genes_for_covariance[gene] = found
        found.append(name)

    def buildRelatedData(self, rsids_from_genes, snps_by_rsid, weights_in_gene):
        related_rsids = []
        related_data = []
        for rsid in rsids_from_genes:
            if not rsid in snps_by_rsid:
                logging.log(5, "related rsid %s not present in genotype data", rsid)
                continue

            related_snp = snps_by_rsid[rsid]
            freq = sum(related_snp.data)*1.0/(2*len(related_snp.data))
            if self.min_maf_filter and self.min_maf_filter > freq:
                logging.log(6, "related rsid %s below min maf: %s", rsid, freq)
                continue

            if self.max_maf_filter and self.max_maf_filter < freq:
                logging.log(6, "related rsid %s  above max maf: %s", rsid, freq)
                continue

            data = related_snp.data
            weight = weights_in_gene[rsid]
            if weight.ref_allele == related_snp.eff_allele and\
                weight.eff_allele == related_snp.ref_allele:
                logging.log(7, "related rsid %s has alleles flipped compared to model, transforming dosage", rsid)
                data = map(lambda x: 2-x, data)

            related_data.append(data)
            related_rsids.append(rsid)
        return related_rsids, related_data

    def buildMatrixOutputEntries(self, matrix, rsids_from_genes, related_rsids, snps_by_rsid):
        entries = []

        #special case: we might have a single rsid!
        if matrix.ndim == 0:
            c = str(float(matrix))
            r = rsids_from_genes[0]
            entries.append((r,r,c))
            return entries

        for i in xrange(0, len(rsids_from_genes)):
            rsid_i = rsids_from_genes[i]
            related_i = -1
            if rsid_i in related_rsids:
                related_i = related_rsids.index(rsid_i)

            for j in xrange(0, len(rsids_from_genes)):
                rsid_j = rsids_from_genes[j]

                related_j = -1
                if rsid_j in related_rsids:
                    related_j = related_rsids.index(rsid_j)

                value = "NA"
                if related_i > -1 and related_j > -1:
                    value = str(matrix[related_i][related_j])

                if i == j:
                    entries.append((rsid_i, rsid_i, value))
                else:
                    if value == "NA":
                        if rsid_i < rsid_j:
                            entries.append((rsid_i, rsid_j, value))
                    else:
                        snp_i = snps_by_rsid[rsid_i]
                        snp_j = snps_by_rsid[rsid_j]

                        if snp_i.position < snp_j.position:
                            entries.append((rsid_i, rsid_j, value))
        return entries

    def writeCovarianceDB(self, output_name, entries):
        connection = sqlite3.connect(output_name)
        cursor = connection.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS covariances (rsid1 CHAR, rsid2 CHAR, covariance TEXT);")
        logging.log(6, "%d covariances inserted in %s", len(entries), output_name)
        cursor.executemany("INSERT INTO covariances(rsid1, rsid2, covariance) VALUES (?,?,?)", entries)
        connection.commit()
        connection.close()

    def addToCorrelationFile(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Building correlation database for %s-%s", name, self.weight_db)
        genes = weight_db_logic.weights_by_gene_name.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            if percent == last_reported_percent+10:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent

            entries = self.buildCorrelationEntries(name, gene, weight_db_logic, snps_by_rsid)

            if len(entries) == 0:
                logging.log(6,"Gene %s has no snps in current file", gene)
                continue

            path = self.correlationFilePath()
            self.addToFile(path, gene, entries)

    def buildCorrelationDBS(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Building correlation database for %s-%s", name, self.weight_db)
        genes = weight_db_logic.weights_by_gene_name.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            if percent == last_reported_percent+10:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent

            output_name_base = "cor-"+name+"-"+gene+"-"+self.weight_db
            output_name = self.correlation_output_folder + "/" + output_name_base
            if os.path.exists(output_name):
                logging.info("Correlation matrix for %s already exists, delete if you want it figured out again", gene)
                continue

            entries = self.buildCorrelationEntries(name, gene, weight_db_logic, snps_by_rsid)

            if len(entries) == 0:
                logging.log(6,"Gene %s has no snps in current file", gene)
                continue

            self.writeCorrelationDB(output_name, entries)

    def buildCorrelationEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        weights_in_gene = weight_db_logic.weights_by_gene_name[gene]
        rsids_from_genes = weights_in_gene.keys()

        #gather as much data as we can work on
        related_rsids, related_data = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)

        if len(related_rsids) == 0:
            return []

        self.updateFoundCorrelation(gene, name)

        #correlation matrix of related SNP's data
        array = numpy.array(related_data)
        cor = numpy.corrcoef(array)

        #translate into sql entries
        entries = self.buildMatrixOutputEntries(cor, rsids_from_genes, related_rsids, snps_by_rsid)
        if not len(entries):
            raise NameError("Couldn not build correlation entries for (%s,%s)" %(name,gene))
        return entries

    def writeCorrelationDB(self, output_name, entries):
        connection = sqlite3.connect(output_name)
        cursor = connection.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS correlations (rsid1 CHAR, rsid2 CHAR, r TEXT);")
        logging.log(6, "%d correlations inserted in %s", len(entries), output_name)
        cursor.executemany("INSERT INTO correlations(rsid1, rsid2, r) VALUES (?,?,?)", entries)
        connection.commit()
        connection.close()

    def updateFoundCorrelation(self, gene, name):
        found = None
        if gene in self.found_genes_for_correlation:
            found = self.found_genes_for_correlation[gene]
            logging.info("Gene %s found again for %s", gene, name)
        else:
            found = []
            self.found_genes_for_correlation[gene] = found
        found.append(name)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build correlations and/or covariances from PHASE3 data and weights database.')

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--weight_db",
                        help="name of weight db in data folder",
                        default="data/DGN-WB_0.5.db")

    parser.add_argument("--input_folder",
                        help="name of folder containing PHASE 3 data",
                        default="intermediate/filtered_1000GP_Phase3")

    parser.add_argument("--correlation_output_folder",
                        help="name of folder to dump results in",
                        default="intermediate/1000GP_Phase3_chr_cor")

    parser.add_argument("--correlations", help="Output correlations",
                        action="store_true",
                        default=False)

    parser.add_argument("--covariance_output_folder",
                        help="name of folder to dump results in",
                        default="intermediate/1000GP_Phase3_chr_cov")

    parser.add_argument("--covariances",
                    help="Output covariances",
                    action="store_true",
                    default=True)

    parser.add_argument('--input_format',
                   help='Input dosage files format. Valid options are: IMPUTE, PrediXcan',
                   default=Formats.PrediXcan)

    parser.add_argument('--output_format',
                   help="Output covariance files format. Valid options are: 'FlatFile', 'DatabaseFile'",
                   default=Formats.FlatFile)

    parser.add_argument('--min_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument('--max_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)


    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = ProcessWeightDB(args)
    work.run()
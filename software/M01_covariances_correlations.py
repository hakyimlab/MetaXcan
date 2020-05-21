#!/usr/bin/env python
__author__ = 'heroico'
import os
import io
import gzip
import logging
import numpy

from metax import WeightDBUtilities
from metax import PrediXcanFormatUtilities
from metax import ThousandGenomesUtilities
from metax import Logging
from metax import Utilities
from metax import Formats
from timeit import default_timer as timer


def pathLeaf(path):
    head, tail = os.path.split(path)
    return tail or os.path.basename(head)

class ProcessWeightDB(object):
    def __init__(self, args):
        self.weight_db = pathLeaf(args.weight_db)
        self.db_path = args.weight_db
        self.data_folder = args.input_folder
        self.correlation_output = args.correlation_output
        self.covariance_output = args.covariance_output
        if args.covariance_output is None:
            comp = os.path.splitext(self.weight_db)[0]
            name = comp + ".cov.txt.gz"
            path = os.path.join("intermediate", "cov")
            path = os.path.join(path, name)
            self.covariance_output = path

        self.input_format = args.input_format

        self.found_genes_for_covariance = {}
        self.found_genes_for_correlation = {}

        self.min_maf_filter = float(args.min_maf_filter) if args.min_maf_filter else None
        self.max_maf_filter = float(args.max_maf_filter) if args.max_maf_filter else None

        self.max_snps_in_gene = int(args.max_snps_in_gene) if args.max_snps_in_gene else None

        self.delimiter = args.delimiter

    def run(self):
        start = timer()

        if not self.correlation_output and not self.covariance_output:
            logging.info("Provide --correlation_output or --covariance_output or both")
            return

        logging.info("Loading Weights")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.db_path)

        logging.info("Building files")
        self.buildFiles(weight_db_logic)


        end = timer()
        logging.info("Ran successfully in %s seconds",  str(end-start))

    def buildFiles(self, weight_db_logic):
        do_correlations = self.correlation_output is not None
        if do_correlations:
            if os.path.exists(self.correlation_output):
                logging.info("%s already exists, delete it if you want it figured out again", self.correlation_output)
                do_correlations = False
            else:
                correlation_dir = os.path.dirname(self.correlation_output)
                if not os.path.exists(correlation_dir):
                    os.makedirs(correlation_dir)
                self.writeFileHeader(self.correlation_output)

        do_covariances = self.covariance_output is not None
        if do_covariances:
            if os.path.exists(self.covariance_output):
                logging.info("%s already exists, delete it if you want it figured out again", self.covariance_output)
                do_covariances = False
            else:
                covariance_dir = os.path.dirname(self.covariance_output)
                if not os.path.exists(covariance_dir):
                    os.makedirs(covariance_dir)
                self.writeFileHeader(self.covariance_output)

        if not do_covariances and not do_correlations:
            return

        names = Utilities.dosageNamesFromFolder(self.data_folder)
        for name in names:
            snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)
            if do_correlations:
                self.addToCorrelationFile(weight_db_logic, name, snps, snps_by_rsid)

            if do_covariances:
                self.addToCovarianceFile(weight_db_logic, name, snps, snps_by_rsid)

    def writeFileHeader(self,path):
        with io.TextIOWrapper(gzip.open(path, "ab"), newline="") as file:
            file.write("GENE RSID1 RSID2 VALUE\n")

    def getSNPS(self, name, weight_db_logic):
        dosageLoader = None
        if self.input_format == Formats.IMPUTE:
            #TODO: update
            dosageLoader = ThousandGenomesUtilities.IMPUTEDosageLoader(self.data_folder, name) #outdated code
        elif self.input_format == Formats.PrediXcan:
            #dosageName = Utilities.dosageName(name)
            path = os.path.join(self.data_folder, name)
            dosageLoader = PrediXcanFormatUtilities.PrediXcanFormatDosageLoader(path, weight_db_logic, self.delimiter)
        else:
            logging.info("Invalid input format: %s", self.input_format)
            return
        snps, snps_by_rsid = dosageLoader.load()
        return snps, snps_by_rsid

    def addToCovarianceFile(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Adding to covariance for %s-%s", name, self.weight_db)

        genes = weight_db_logic.weights_by_gene.keys()
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

            self.addToFile(self.covariance_output, gene, entries)

    def addToFile(self, path, gene, entries):
        with io.TextIOWrapper(gzip.open(path, "ab"), newline="") as file:
            for entry in entries:
                line = " ".join([gene, entry[0], entry[1], entry[2]])+"\n"
                file.write(line)

    def buildCovarianceEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = list(weights_in_gene.keys())

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

        l = len(rsids_from_genes)
        if self.max_snps_in_gene and l > self.max_snps_in_gene:
            logging.info("Skipping covariance too large: %d", l)
            return related_data, related_rsids

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
                data = [2-x for x in data]

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

        for i in range(0, len(rsids_from_genes)):
            rsid_i = rsids_from_genes[i]
            related_i = -1
            if rsid_i in related_rsids:
                related_i = related_rsids.index(rsid_i)

            for j in range(0, len(rsids_from_genes)):
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

    def addToCorrelationFile(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Building correlation database for %s-%s", name, self.weight_db)
        genes = list(weight_db_logic.weights_by_gene.keys())
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

            self.addToFile(self.correlation_output, gene, entries)

    def buildCorrelationEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = list(weights_in_gene.keys())

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
                        default="intermediate/TGF_EUR")

    parser.add_argument("--correlation_output",
                        help="Name of file to dump correlation results in.",
                        default=None)
                        #default="intermediate/1000GP_Phase3_chr_cor_dgnwb_cor.txt.gz")

    parser.add_argument("--covariance_output",
                        help="Name of file to dump covariance results in. Defaults to 'intermediate/cov/' + file name prefix from '--weight_db' argument",
                        default=None)

    parser.add_argument('--input_format',
                   help='Input dosage files format. Valid options are: IMPUTE, PrediXcan',
                   default=Formats.PrediXcan)

    parser.add_argument('--min_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument('--max_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument("--max_snps_in_gene",
                        help="Ignore any gene that has snps above this value",
                        type=int,
                        default=None)

    parser.add_argument("--delimiter", help="dosage delimiter", default=" ")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = ProcessWeightDB(args)
    work.run()
#! /usr/bin/env python
__author__ = 'heroico'

import logging
import numpy
import gzip
import os
import metax.WeightDBUtilities as WeightDBUtilities
import metax.Utilities as Utilities
import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.Logging as Logging

class CalculateVariances(object):
    def __init__(self, args):
        self.weight_db = args.weight_db
        self.data_folder_phase = args.phase_folder
        self.output_file = args.output_file

    def run(self):
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.weight_db) if self.weight_db else None
        contents = Utilities.contentsWithPatternsFromFolder(self.data_folder_phase, ["gz"])

        if os.path.exists(self.output_file):
            logging.info("Variance output already exists, delete it if you want stuff to be figured out again")
            return

        dir = os.path.dirname(self.output_file)
        if not os.path.exists(dir):
            os.makedirs(dir)

        for content in contents:
            self.buildVarianceDB(weight_db_logic,content)

    def buildVarianceDB(self, weight_db_logic, content):
        logging.info("Building variance database for %s-%s", content, (self.weight_db if self.weight_db else "No db"))
        path = os.path.join(self.data_folder_phase, content)
        dosageLoader = PrediXcanFormatUtilities.PrediXcanFormatDosageLoader(path, weight_db_logic)
        snps, snps_by_rsid = dosageLoader.load()

        total_snps = len(snps)
        last_reported_percent = 1
        processed = 0

        with gzip.open(self.output_file, 'ab') as file:
            for snp in snps:
                processed += 1
                percent = processed*100.0 / total_snps
                if percent > last_reported_percent+9:
                    logging.info("%d percent snps processed", percent)
                last_reported_percent = percent

                if weight_db_logic and not snp.name in weight_db_logic.genes_for_an_rsid:
                    logging.log(6, "rsid %s not in tissue database", snp.name)
                    continue

                var = numpy.var(snp.data, ddof=1)

                line = ",".join([snp.name, str(float(var))])+"\n"
                file.write(line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build variances from PHASE3 data and weights database.')

    parser.add_argument("--weight_db",
                        help="name of weight db in data folder")
                        #,default="data_folder/DGN-WB_0.5.db")

    parser.add_argument("--phase_folder",
                        help="name of folder containing PHASE 3 data",
                        default="intermediate/TGF_EUR")

    parser.add_argument("--output_file",
                        help="name of folder to dump results in",
                        default="intermediate/var/var.txt.gz")

    args = parser.parse_args()

    Logging.configureLogging(7)

    work = CalculateVariances(args)
    work.run()


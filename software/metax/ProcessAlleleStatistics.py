#!/usr/bin/env python
__author__ = 'heroico'

import logging
import os
import ThousandGenomesUtilities
import GWASUtilities
import Utilities
import WeightDBUtilities
import Logging


class AlleleStats(object):
    def __init__(self, rsid, chromosome, weight_db_entry):
        self.rsid = rsid
        self.chromosome = chromosome

        self.weight_db_ref_allele = weight_db_entry.ref_allele
        self.weight_db_eff_allele = weight_db_entry.eff_allele

        self.gwas_ref_allele = "NA"
        self.gwas_eff_allele = "NA"
        self.gwas_OR_BETA = "NA"

        self.legend_type = "NA"
        self.legend_ref_allele = "NA"
        self.legend_eff_allele = "NA"

    def getGWASEntryData(self, gwas_entry):
        if gwas_entry:
            self.gwas_ref_allele = gwas_entry[GWASUtilities.GWASSNPInfoLineCollector.A1]
            self.gwas_eff_allele = gwas_entry[GWASUtilities.GWASSNPInfoLineCollector.A2]
            self.gwas_OR_BETA = gwas_entry[GWASUtilities.GWASSNPInfoLineCollector.OR_BETA]
        else:
            self.gwas_ref_allele = "NA"
            self.gwas_eff_allele = "NA"
            self.gwas_OR_BETA = "NA"

    def getLegendData(self, type, a0, a1):
        self.legend_type = type if type is not None else "NA"
        self.legend_ref_allele = a0 if a0 is not None else "NA"
        self.legend_eff_allele = a1 if a1 is not  None else "NA"

    def toCSVLine(self):
        tuple = (self.rsid, self.chromosome,
                 self.weight_db_ref_allele, self.weight_db_eff_allele,
                 self.legend_ref_allele, self.legend_eff_allele, self.legend_type,
                 self.gwas_ref_allele, self.gwas_eff_allele, self.gwas_OR_BETA)
        line = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % tuple
        return line

    @classmethod
    def CSVHeader(cls):
        return "rsid,chromosome,wdb_ref_allele,wdb_eff_allele,legend_ref_allele,legend_eff_allele,legend_type,gwas_ref_allele,gwas_eff_allele,gwas_OR_BETA\n"

class ProcessAlleleStatistics(object):
    def __init__(self, args):
        self.data_folder = args.data_folder

        self.weight_db = args.weight_db
        self.db_path = os.path.join(self.data_folder, args.weight_db)

        self.data_folder_phase = args.phase_folder
        self.data_folder_gwas_dosage = args.gwas_dosage_folder

        self.output_file = args.output_file

    def run(self):
        if os.path.exists(self.output_file):
            logging.info("File %s already exists, delete it if you want it calculated again", self.output_file)
            return

        logging.info("Opening %s", self.weight_db)
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.db_path)

        CHROMOSOMES = ["chr"+str(x) for x in xrange(1, 23)]

        dosage_names = Utilities.dosageNamesFromFolder(self.data_folder_gwas_dosage)
        legend_names = Utilities.legendNamesFromFolder(self.data_folder_phase)

        findings={}
        for chromosome in CHROMOSOMES:
            logging.info("Processing chromosome %s", chromosome)
            dosage_name = Utilities.removeNameWithPatterns(dosage_names, [chromosome+"."])
            dosage = self.loadDosageFile(self.data_folder_gwas_dosage, dosage_name)
            self.processDosage(chromosome, weight_db_logic, dosage, findings)

            legend_name = Utilities.removeNameEndingWith(legend_names, chromosome)
            self.processLegendName(chromosome, weight_db_logic, dosage, findings, legend_name)


        with open(self.output_file, "w") as file:
            file.write(AlleleStats.CSVHeader())

            def sortByChromosome(finding):
                return finding.chromosome
            entries = sorted(findings.values(), key=sortByChromosome)
            for finding in entries:
                line = finding.toCSVLine()
                file.write(line)

    def loadDosageFile(self, base_path, name):
        callback = GWASUtilities.GWASSNPInfoLineCollector()
        dosage_loader = GWASUtilities.GWASDosageFileLoader(base_path, name,  callback)
        keyed_data_set = dosage_loader.load()
        return keyed_data_set

    def processDosage(self, chromosome, weight_db_logic, dosage, findings):
        ok = 0
        for rsid, dosage_entry in dosage.values_by_key.iteritems():
            weight_db_entry = weight_db_logic.anEntryWithRSID(rsid)
            if not weight_db_entry:
                logging.log(7, "%s in dosage not in weights", rsid)
                continue

            a1 = dosage_entry[GWASUtilities.GWASSNPInfoLineCollector.A1]
            a2 = dosage_entry[GWASUtilities.GWASSNPInfoLineCollector.A2]
            OR= dosage_entry[GWASUtilities.GWASSNPInfoLineCollector.OR]
            if not weight_db_entry.ref_allele == a1 or \
                not weight_db_entry.eff_allele == a2 or \
                OR == "NA":
                logging.log(7, "%s in dosage is problematic (%s, %s)(%s, %s, %s)", rsid, weight_db_entry.ref_allele, weight_db_entry.eff_allele, a1, a2, OR)
                finding = appropriateFinding(findings, rsid, chromosome, weight_db_entry)
                finding.getGWASEntryData(dosage_entry)
                continue

            ok += 1
        logging.log(8,"After processing dosage, %d snps were found to be ok", ok)


    def processLegendName(self, chromosome, weight_db_logic, dosage, findings, legend_name):
        class LegendCallback(object):
            def __init__(self):
                pass

            def __call__(self, i, comps):
                id = comps[ThousandGenomesUtilities.ILTF.ID]
                id_components = id.split(':')
                rsid = id_components[0]

                weight_db_entry = weight_db_logic.anEntryWithRSID(rsid)
                if not weight_db_entry:
                    logging.log(8, "rsid %s from legend not in db, %s", rsid, id)
                    return

                type = comps[ThousandGenomesUtilities.ILTF.TYPE]
                a0 = comps[ThousandGenomesUtilities.ILTF.A0]
                a1 = comps[ThousandGenomesUtilities.ILTF.A1]

                if rsid in findings:
                    finding = findings[rsid]
                    finding.getLegendData(type, a0, a1)

                move_on = True
                if not type == "Biallelic_SNP":
                    logging.log(8, "%s %s Not biallelic: %s", chromosome, id, type)
                    move_on = False
                else:
                    if (a0 == 'T' and a1 == 'A') or \
                        (a0 == 'A' and a1 == 'T') or \
                        (a0 == 'C' and a1 == 'G') or \
                        (a0 == 'G' and a1 == 'C'):
                        logging.log(8, "%s %s Problematic: %s, %s", chromosome, id, a0, a1)
                        move_on = False

                    if not weight_db_entry.ref_allele == a0 or \
                        not weight_db_entry.eff_allele == a1:
                        logging.log(8, "%s %s Different alleles %s %s", chromosome, id, a0, a1)
                        move_on = False

                if not move_on:
                    finding = appropriateFinding(findings, rsid, chromosome, weight_db_entry)
                    finding.getLegendData(type,a0,a1)

                    dosage_entry = None
                    if rsid in dosage.values_by_key:
                        dosage_entry = dosage.values_by_key[rsid]
                    finding.getGWASEntryData(dosage_entry)



        callback = LegendCallback()
        loader = ThousandGenomesUtilities.LEGENDLoader(self.data_folder_phase, legend_name)
        loader.iterateOverFileLegends(callback)

def appropriateFinding(findings, rsid, chromosome, weight_db_entry):
    finding = None
    if rsid in findings:
        finding = findings[rsid]
    else:
        finding = AlleleStats(rsid, chromosome, weight_db_entry)
        findings[rsid] = finding
    return finding

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build correlations from PHASE3 data and weights database.')

    parser.add_argument("--data_folder",
                        help="higher level data folder",
                        default="data")

    parser.add_argument("--weight_db",
                        help="name of weight db in data folder",
                        default="DGN-WB_0.5.db")

    parser.add_argument("--phase_folder",
                        help="name of folder containing PHASE 3 data",
                        default="data/1000GP_Phase3")

    parser.add_argument("--gwas_dosage_folder",
                        help="name of folder containing dosage data",
                        default="data/T1D-GWAS")

    parser.add_argument("--output_file",
                        help="name of file to dump results in",
                        default="results/allele_stats.csv")



    args = parser.parse_args()

    Logging.configureLogging(logging.INFO)

    work = ProcessAlleleStatistics(args)
    work.run()

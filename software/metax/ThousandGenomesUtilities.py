__author__ = 'heroico'

import gzip
import os
import io
import logging
from . import DataSetSNP
from . import Utilities

class ILTF:
    """IMPUTE legend file format"""
    ID = 0
    POSITION = 1
    A0 = 2
    A1 = 3
    TYPE = 4
    AFR = 5
    AMR = 6
    EAS = 7
    EUR = 8
    SAS = 9
    ALL = 10

class LEGENDLoader(object):
    def __init__(self, base_path=None, name=None):
        self.base_path = base_path
        self.name = name
        self.chrlegend = Utilities.legendName(name)
        self.chrnumber = name.split('chr')[1]

    def iterateOverFileLegends(self, callback):
        path = os.path.join(self.base_path, self.chrlegend)
        file_iterator = Utilities.CSVFileIterator(path, header="id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL", compressed=True)
        file_iterator.iterate(callback)

class IMPUTELoader(object):
    def __init__(self, base_path=None, name=None):
        self.base_path = base_path
        self.name = name
        self.chrhap = Utilities.hapName(name)
        self.chrlegend = Utilities.legendName(name)
        self.chrdosage = Utilities.dosageName(name)
        self.chrnumber = name.split('chr')[1]

    def iterateOverFile(self, row_processor):
        """row_processor should provide an appropiate '__call__' method"""
        with io.TextIOWrapper(gzip.open(os.path.join(self.base_path, self.chrlegend), 'r'), newline="") as legend_file:
            legend_file.readline()
            with io.TextIOWrapper(gzip.open(os.path.join(self.base_path, self.chrhap), 'r'), newline="") as hap_file:
                row=1
                last_reported_row = 1
                while True:
                    logging.log(2,"Processing row %d", row)
                    row += 1
                    if row - last_reported_row > 999999:
                        logging.log(8, "Processed row %d", row)
                        last_reported_row = row

                    legend_line = legend_file.readline()
                    hap_line = hap_file.readline()
                    if hap_line == '' or legend_line == None:
                        break

                    row_processor(hap_line, legend_line)

    def iterateOverFileDosage(self, row_processor):
        """row_processor should provide an appropiate '__call__' method"""
        with io.TextIOWrapper(gzip.open(os.path.join(self.base_path, self.chrlegend), 'r'), newline="") as legend_file:
            legend_file.readline()
            with io.TextIOWrapper(gzip.open(os.path.join(self.base_path, self.chrdosage), 'r'), newline="") as dosage_file:
                row=1
                last_reported_row = 1
                while True:
                    logging.log(2,"Processing row %d", row)

                    legend_line = legend_file.readline()
                    dosage_line = dosage_file.readline()
                    if legend_line == '' or legend_line == None or dosage_line == '' or dosage_line == None:
                        break

                    row_processor(dosage_line, legend_line, row-1)

                    row += 1
                    if row - last_reported_row > 99999:
                        logging.log(8, "Processed row %d", row)
                        last_reported_row = row

class IMPUTEDosageLoader(object):
    def __init__(self, base_path = None, name=None):
        self.base_path = base_path
        self.name = name

    def load(self):
        logging.info("Loading %s dosage from %s", self.name, self.base_path)
        class SNPCollector(object):
            def __init__(self, snps=[], snps_by_rsid={}):
                self.snps = snps
                self.snps_by_rsid = snps_by_rsid

            def __call__(self, dosage_line, legend_line, row):
                legend = legend_line.split(" ")
                id = legend[ILTF.ID]
                id_components = id.split(':')
                rsid = id_components[0]
                position = id_components[1]
                ref_allele = legend_line[ILTF.A0]
                eff_allele = legend_line[ILTF.A1]
                data = list(map(int, dosage_line.strip().split(" ")))

                snp = DataSetSNP.DataSetSNP(name=rsid, index=row, data=data, position=int(position), ref_allele=ref_allele, eff_allele=eff_allele)
                if snp.name in self.snps_by_rsid:
                    old = self.snps_by_rsid[snp.name]
                    logging.info("Duplicated rsid: (%s,%s) %s", old.name, old.position, legend_line.strip())

                self.snps.append(snp)
                self.snps_by_rsid[snp.name] = snp

        loader = IMPUTELoader(self.base_path,self.name)
        snp_collector = SNPCollector()
        loader.iterateOverFileDosage(snp_collector)
        return snp_collector.snps, snp_collector.snps_by_rsid

class IMPUTEFilteredDosageFileBuilder(object):
    """Builds a filtered file selecting from rsid and people"""
    def __init__(self, base_path = None, name = None, output_pattern = None, snp_dict = None, all_people=None, selected_people_by_id=None, chromosome_name=None):
        self.base_path = base_path
        self.name = name
        self.output_pattern = output_pattern
        self.chromosome_name = chromosome_name
        self.all_people = all_people if all_people is not None else []
        self.selected_people_by_id = selected_people_by_id if selected_people_by_id is not None else {}
        self.snp_dict = snp_dict if snp_dict is not None else {}

    def buildIMPUTE(self):
        logging.info("Building IMPUTE dosage file for %s", os.path.join(self.base_path,self.name))
        logging.debug("all people: %d", len(self.all_people))
        logging.debug("selected people: %d", len(list(self.selected_people_by_id.keys())))
        logging.debug("snps: %d", len(list(self.snp_dict.keys())))

        dosage_file_name = Utilities.dosageName(self.output_pattern)
        legend_file_name = Utilities.legendName(self.output_pattern)
        if os.path.exists(dosage_file_name) and os.path.exists(legend_file_name):
            logging.info("Files for %s already built, delete it if you want it done again", dosage_file_name)
            return

        class IMPUTEOutput(object):
            def __init__(self, snp_dict, all_people, selected_people_by_id, hap_output_file, legend_output_file):
                self.hap_output_file = hap_output_file
                self.legend_output_file = legend_output_file
                self.snp_dict = snp_dict
                self.all_people = all_people
                self.selected_people_by_id = selected_people_by_id
                self.found_snps = 0
                self.last_reported_percenteage = 0
                self.snps_len = len(list(self.snp_dict.keys()))

            def __call__(self, hap_line, legend_line):
                rsid, valid, legend = checkLegend(legend_line, self.snp_dict)
                if not valid:
                    return

                logging.log(3, "id %s found in white list", rsid)
                self.found_snps += 1
                percent = int(round(self.found_snps * 100.0 / self.snps_len))
                if percent > self.last_reported_percenteage:
                    self.last_reported_percenteage = percent
                    logging.log(9, "%d percent of snps found", percent)

                self.legend_output_file.write(legend_line)

                dosages = buildDosages(hap_line, self.all_people, self.selected_people_by_id)
                dosages_line = " ".join(map(str,dosages))+"\n"
                self.hap_output_file.write(dosages_line)

        with gzip.open(dosage_file_name, 'wb') as hap_output_file:
            with gzip.open(legend_file_name, 'wb') as legend_output_file:
                legend_output_file.write("id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL\n")
                loader = IMPUTELoader(self.base_path, self.name)
                callback = IMPUTEOutput(self.snp_dict, self.all_people, self.selected_people_by_id, hap_output_file, legend_output_file)
                loader.iterateOverFile(callback)

    def buildPrediXcan(self):
        logging.info("Building PrediXcan dosage file for %s", os.path.join(self.base_path,self.name))
        logging.debug("all people: %d", len(self.all_people))
        logging.debug("selected people: %d", len(list(self.selected_people_by_id.keys())))
        logging.debug("snps: %d", len(list(self.snp_dict.keys())))

        dosage_file_name = Utilities.dosageName(self.output_pattern)
        if os.path.exists(dosage_file_name):
            logging.info("Files for %s already built, delete it if you want it done again", dosage_file_name)
            return

        class PrediXcanOutput(object):
            def __init__(self, snp_dict, all_people, selected_people_by_id, dosage_output_file, chromosome_name):
                self.chromosome_name = chromosome_name
                self.dosage_output_file = dosage_output_file
                self.snp_dict = snp_dict
                self.all_people = all_people
                self.selected_people_by_id = selected_people_by_id
                self.found_snps = 0
                self.last_reported_percenteage = 0
                self.snps_len = len(list(self.snp_dict.keys()))

            def __call__(self, hap_line, legend_line):
                rsid, valid, legend = checkLegend(legend_line, self.snp_dict)
                if not valid:
                    return

                logging.log(3, "id %s found in white list", rsid)
                self.found_snps += 1
                percent = int(round(self.found_snps * 100.0 / self.snps_len))
                if percent > self.last_reported_percenteage:
                    self.last_reported_percenteage = percent
                    logging.log(9, "%d percent of snps found", percent)

                dosages = buildDosages(hap_line, self.all_people, self.selected_people_by_id)
                snp_string = buildPrediXcanSNPFields(legend, dosages, self.chromosome_name)
                dosages_line = " ".join(map(str,dosages))+"\n"
                line = snp_string + " " + dosages_line
                self.dosage_output_file.write(line)

        with io.TextIOWrapper(gzip.open(dosage_file_name, 'w'), newline="") as dosage_output_file:
            loader = IMPUTELoader(self.base_path, self.name)
            callback = PrediXcanOutput(self.snp_dict, self.all_people, self.selected_people_by_id, dosage_output_file, self.chromosome_name)
            loader.iterateOverFile(callback)

def checkLegend(legend_line, snp_dict):
    legend = legend_line.split(" ")
    id = legend[ILTF.ID]
    rsid = id.split(':')[0]
    logging.log(3, "%s",rsid)
    if not rsid in snp_dict:
        logging.log(3, "id not found in white list, skipping")
        return rsid, False, legend

    a0 = legend[ILTF.A0]
    a1 = legend[ILTF.A1]
    if (a0 == 'T' and a1 == 'A') or \
        (a0 == 'A' and a1 == 'T') or \
        (a0 == 'C' and a1 == 'G') or \
        (a0 == 'G' and a1 == 'C'):
        logging.log(4, "discarding problematic case")
        return rsid, False, legend

    type = legend[ILTF.TYPE]
    if type != 'Biallelic_SNP':
        logging.log(4, "Discarding rsid with dubious type")
        return rsid, False, legend

    return rsid, True, legend

def buildPrediXcanSNPFields(legend, dosage, chromosome_name):
    average = float(sum(dosage))/(2*len(dosage))
    average = str(average)

    id = legend[ILTF.ID]
    id_comps = id.split(':')
    rsid = id_comps[0]
    position = id_comps[1]
    a0 = legend[ILTF.A0]
    a1 = legend[ILTF.A1]

    fields = " ".join([chromosome_name, rsid, position, a0, a1, average])
    return fields

def buildDosages(hap_line, all_people, selected_people_by_id):
    dosages = []
    hap = hap_line.split(" ")
    len_hap = len(hap)
    len_all = len(all_people)
    if len_hap != 2*len_all:
        logging.error("not enough entries: (%d,%d)", len_hap, len_all)
        assert False

    for index,person in enumerate(all_people):
        if not person.id in selected_people_by_id:
            continue
        a0 = int(hap[2*index])
        a1 = int(hap[2*index+1])
        dosages.append(a0+a1)
    return dosages

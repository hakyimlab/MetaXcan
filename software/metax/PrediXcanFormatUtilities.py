__author__ = 'heroico'

import logging
import os
import gzip
from . import Utilities
from . import DataSetSNP


class PDTF:
    """Format of PrediXcan dosage"""
    CHR = 0
    RSID = 1
    POSITION = 2
    ALLELE_0 = 3
    ALLELE_1 = 4
    COLUMN_5 = 5
    FIRST_DATA_COLUMN = 6

class PrediXcanFormatDosageLoader(object):
    #weight_db_logic is used only for discarding absent snps.
    def __init__(self, path, weight_db_logic, delimiter = " "):
        self.path = path
        self.weight_db_logic = weight_db_logic
        self.delimiter = delimiter

    def load(self):
        logging.info("Loading %s dosage", self.path)
        class PrediXcanCollector(object):
            def __init__(self, snps=[], snps_by_rsid={}, weight_db_logic=None):
                self.snps = snps
                self.snps_by_rsid = snps_by_rsid
                self.weight_db_logic = weight_db_logic

            def __call__(self, i, components):
                rsid = components[PDTF.RSID]
                if self.weight_db_logic and not rsid in self.weight_db_logic.genes_for_an_rsid:
                    logging.log(5, "rsid %s not in weight db, skip it", rsid)
                    return

                position = components[PDTF.POSITION]

                ref_allele = components[PDTF.ALLELE_0]
                if not ref_allele in Utilities.VALID_ALLELES:
                    logging.log(9, "wrong ref allele, rsid %s is not an SNP", rsid)
                    return
                eff_allele = components[PDTF.ALLELE_1]
                if not eff_allele in Utilities.VALID_ALLELES:
                    logging.log(9, "wrong eff allele, rsid %s is not an SNP", rsid)
                    return
                dosages = list(map(float,components[PDTF.FIRST_DATA_COLUMN:])) #dosages may be inputed
                #Should we flip based on weight_db at this point?

                snp = DataSetSNP.DataSetSNP(name=rsid, index=i, data=dosages, position=int(position), ref_allele=ref_allele, eff_allele=eff_allele)
                if snp.name in self.snps_by_rsid:
                    old = self.snps_by_rsid[snp.name]
                    logging.log(9, "Duplicated rsid: (%s,%s)", old.name, old.position)
                    return
                self.snps.append(snp)
                self.snps_by_rsid[snp.name] = snp

        loader = Utilities.CSVFileIterator(self.path, compressed=True, delimiter=self.delimiter)
        collector = PrediXcanCollector(weight_db_logic=self.weight_db_logic)
        loader.iterate(collector)
        return collector.snps, collector.snps_by_rsid

class PrediXcanFormatFilteredFilesProcess(object):
    def __init__(self, input_file_path, output_path, name, all_people, selected_people_by_id, snps_dict, delimiter=" "):
        self.input_file_path = input_file_path
        self.output_path = output_path
        self.name = name
        self.all_people = all_people
        self.selected_people_by_id = selected_people_by_id
        self.snps_dict = snps_dict
        self.delimiter = delimiter

    def buildIMPUTE(self):
        dosage_path = dosagePath(self.output_path, self.name)
        legend_path = legendPath(self.output_path, self.name)
        if os.path.exists(dosage_path) or os.path.exists(legend_path):
            logging.info("%s and/or %s already exists, delete it if you want it done again", dosage_path, legend_path)
            return

        class IMPUTEOutput(object):
            def __init__(self, snps_dict, all_people, selected_people_by_id, legend_file, dosage_file):
                self.all_people = all_people
                self.selected_people_by_id = selected_people_by_id
                self.snps_dict = snps_dict
                self.legend_file = legend_file
                self.dosage_file = dosage_file

            def __call__(self, i, row):
                rsid = row[PDTF.RSID]

                if rsid not in self.snps_dict:
                    logging.log(5, "rsid %s not in whitelist", rsid)
                    return

                position = row[PDTF.POSITION]
                a0 = row[PDTF.ALLELE_0]
                a1 = row[PDTF.ALLELE_1]

                first = "%s:%s:%s:%s" % (rsid, position, a0, a1)
                fields = [first, position, a0, a1, "Biallelic_SNP", "NA", "NA", "NA", "NA", "NA", "NA"]
                legend_line = " ".join(fields)+"\n"
                self.legend_file.write(legend_line)


                dosages = row[PDTF.FIRST_DATA_COLUMN:]
                if len(dosages) != len(self.all_people):
                    logging.error("rsid %s: not enough dosage: %d, %d", rsid, len(dosages), len(self.all_people))
                    assert False
                selected_dosages = pickDosages(dosages, self.all_people, self.selected_people_by_id)
                dosage_line = " ".join(selected_dosages)+"\n"
                self.dosage_file.write(dosage_line)

        logging.info("building %s", dosage_path)
        with gzip.open(legend_path, 'wb') as legend_file:
            legend_file.write("id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL\n")
            with gzip.open(dosage_path, 'wb') as dosage_file:
                iterator = Utilities.CSVFileIterator(self.input_file_path, compressed=True)
                callback = IMPUTEOutput(self.snps_dict, self.all_people, self.selected_people_by_id, legend_file, dosage_file)
                iterator.iterate(callback)


    def buildPrediXcan(self):
        dosage_path = dosagePath(self.output_path, self.name)
        if os.path.exists(dosage_path):
            logging.info("%s already exists, delete it if you want it done again", dosage_path)
            return

        class PrediXcanOutput(object):
            def __init__(self, snps_dict, all_people, selected_people_by_id, dosage_file):
                self.all_people = all_people
                self.selected_people_by_id = selected_people_by_id
                self.snps_dict = snps_dict
                self.dosage_file = dosage_file

            def __call__(self, i, row):
                rsid = row[PDTF.RSID]

                if rsid not in self.snps_dict:
                    logging.log(5, "rsid %s not in whitelist", rsid)
                    return

                chromosome = row[PDTF.CHR]
                position = row[PDTF.POSITION]
                a0 = row[PDTF.ALLELE_0]
                a1 = row[PDTF.ALLELE_1]

                dosages = row[PDTF.FIRST_DATA_COLUMN:]
                if len(dosages) != len(self.all_people):
                    logging.error("rsid %s: not enough dosage: %d, %d", rsid, len(dosages), len(self.all_people))
                    assert False
                selected_dosages = pickDosages(dosages, self.all_people, self.selected_people_by_id)
                dosages = list(map(float,selected_dosages)) # dosages may be inputed
                average = float(sum(dosages))/(2*len(dosages))
                average = str(average)

                fields = " ".join([chromosome, rsid, position, a0, a1, average])
                dosage_line = " ".join(selected_dosages)+"\n"
                line = fields + " " + dosage_line
                self.dosage_file.write(line)

        logging.info("building %s", dosage_path)
        with gzip.open(dosage_path, 'wb') as dosage_file:
            iterator = Utilities.CSVFileIterator(self.input_file_path, compressed=True)
            callback = PrediXcanOutput(self.snps_dict, self.all_people, self.selected_people_by_id, dosage_file)
            iterator.iterate(callback)

def pickDosages(dosages, all_people, selected_people_by_id):
    selected = []
    for i, person in enumerate(all_people):
        if not person.id in selected_people_by_id:
            continue
        selected.append(dosages[i])
    return selected

def dosagePath(output_path, name):
    base_name = name.strip(".dosage.txt.gz")
    dosage_name = Utilities.dosageName(base_name)
    output_path = os.path.join(output_path, dosage_name)
    return output_path

def legendPath(output_path, name):
    base_name = name.strip(".dosage.txt.gz")
    legend_name = Utilities.legendName(base_name)
    output_path = os.path.join(output_path, legend_name)
    return output_path

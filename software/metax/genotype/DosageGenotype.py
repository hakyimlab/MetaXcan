import gzip
import os
import logging

import numpy

from Genotype import GF
from .. import Utilities

class DTF:
    """Format of dosage"""
    CHR = 0
    RSID = 1
    POSITION = 2
    ALLELE_0 = 3
    ALLELE_1 = 4
    FREQ= 5
    FIRST_DATA_COLUMN = 6

def dosage_file_geno_lines(file, snps=None):
    logging.log(9, "Processing Dosage geno %s", file)
    with gzip.open(file) as file:
        for i,line in enumerate(file):
            comps = line.strip().split()
            rsid = comps[DTF.RSID]

            if snps and not rsid in snps:
                continue

            dosage = numpy.array(comps[DTF.FIRST_DATA_COLUMN:], dtype=numpy.float64)
            chrom = comps[DTF.CHR].replace("chr", "")

            yield [rsid, chrom, comps[DTF.POSITION], comps[DTF.ALLELE_0], comps[DTF.ALLELE_1], comps[DTF.FREQ]] + list(dosage)

def dosage_folder_geno_lines(dosage_folder, dosage_pattern, snps=None):
    logging.log(9, "Processing Dosage Folder geno")
    files = Utilities.contentsWithRegexpFromFolder(dosage_folder, dosage_pattern)
    files = sorted([os.path.join(dosage_folder, x) for x in files])
    for f in files:
        for line in dosage_file_geno_lines(f, snps=snps):
            yield line

def dosage_geno_by_chromosome(dosage_folder, dosage_pattern, snps=None):
    buffer = []
    last_chr = None

    def _buffer_to_data(buffer):
        _data = zip(*buffer)

        _metadata = _data[0:GF.FIRST_DOSAGE]
        rsids = _metadata[GF.RSID]
        chromosome = int(_metadata[GF.CHROMOSOME][0]) # I know I am gonna regret this cast
        _metadata = zip(*_metadata)
        metadata = Utilities.to_dataframe(_metadata, ["rsid", "chromosome", "position", "ref_allele", "alt_allele", "frequency"], to_numeric="ignore")

        dosage = zip(*_data[GF.FIRST_DOSAGE:])
        dosage_data = {rsids[i]:dosage[i] for i in xrange(0, len(rsids))}

        return chromosome, metadata, dosage_data

    logging.log(8, "Starting to process lines")
    for line in dosage_folder_geno_lines(dosage_folder, dosage_pattern, snps):
        chromosome = line[GF.CHROMOSOME]
        if last_chr is None: last_chr = chromosome

        if last_chr != chromosome:
            yield _buffer_to_data(buffer)
            buffer = []

        last_chr = chromosome
        buffer.append(line)

    yield _buffer_to_data(buffer)


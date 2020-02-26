import re
import os
import logging

import numpy

from .Genotype import GF
from .. import Utilities
from ..misc import Genomics

class DTF:
    """Format of dosage"""
    CHR = 0
    ID = 1
    POSITION = 2
    ALLELE_0 = 3
    ALLELE_1 = 4
    FREQ= 5
    FIRST_DATA_COLUMN = 6

def dosage_file_geno_lines(file, snps=None, skip_palindromic=False):
    logging.log(9, "Processing Dosage geno %s", file)
    for i, line in enumerate(Utilities.generate_from_any_plain_file(file)):
        comps = line.strip().split()
        id = comps[DTF.ID]

        if snps and not id in snps:
            continue

        ref_allele, alt_allele = comps[DTF.ALLELE_0], comps[DTF.ALLELE_1]
        if skip_palindromic and Genomics.is_palindromic(ref_allele, alt_allele):
            continue

        dosage = numpy.array(comps[DTF.FIRST_DATA_COLUMN:], dtype=numpy.float64)
        chrom = comps[DTF.CHR].replace("chr", "")

        yield (id, int(chrom), int(comps[DTF.POSITION]), ref_allele, alt_allele, float(comps[DTF.FREQ])) + tuple(dosage)

def dosage_files_geno_lines(dosage_files, snps=None, skip_palindromic=False):
    chr_ = re.compile(".*chr(\d+).*")
    def sort_geno(x):
        x_ = chr_.search(x).group(1)
        if x_:
            return int(x_)
        else:
            return x
    dosage_files = sorted(dosage_files, key=sort_geno)

    for f in dosage_files:
        for e in dosage_file_geno_lines(f, snps=snps, skip_palindromic=skip_palindromic):
            yield e

def dosage_folder_geno_lines(dosage_folder, dosage_pattern, snps=None):
    logging.log(9, "Processing Dosage Folder geno %s", dosage_folder)
    files = Utilities.contentsWithRegexpFromFolder(dosage_folder, dosage_pattern)
    files = sorted([os.path.join(dosage_folder, x) for x in files])
    for e in dosage_files_geno_lines(files, snps=snps):
        yield e

def dosage_geno_by_chromosome(dosage_folder, dosage_pattern, snps=None):
    buffer = []
    last_chr = None

    def _buffer_to_data(buffer):
        _data = list(zip(*buffer))

        _metadata = _data[0:GF.FIRST_DOSAGE]
        rsids = _metadata[GF.RSID]
        chromosome = int(_metadata[GF.CHROMOSOME][0]) # I know I am gonna regret this cast
        _metadata = list(zip(*_metadata))
        metadata = Utilities.to_dataframe(_metadata, ["rsid", "chromosome", "position", "ref_allele", "alt_allele", "frequency"], to_numeric="ignore")

        dosage = list(zip(*_data[GF.FIRST_DOSAGE:]))
        dosage_data = {rsids[i]:dosage[i] for i in range(0, len(rsids))}

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


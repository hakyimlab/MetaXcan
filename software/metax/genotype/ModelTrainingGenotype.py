__author__ = "alvaro barbeira"

import gzip
import os
import logging

import numpy

from .Genotype import GF
from .. import Utilities
from . import Helpers

def model_training_file_geno_lines(file, snp_annotation, snps=None):
    logging.log(9, "Processing Dosage geno %s", file)
    with gzip.open(file) as file:
        for i,line in enumerate(file):
            if i==0: continue
            comps = line.strip().split()
            id = comps[0]

            if not id  in snp_annotation:
                continue

            rsid = snp_annotation[id]
            if snps and not rsid in snps:
                continue

            dosage = numpy.array(comps[1:], dtype=numpy.float64)
            chrom, position, allele_0, allele_1, leftover = id.split("_")

            yield [rsid, chrom, position, allele_0, allele_1, numpy.mean(dosage)/2] + list(dosage)

def model_training_folder_geno_lines(dosage_folder, dosage_pattern, snp_annotation, snps=None):
    logging.log(9, "Processing Dosage Folder geno")
    files = Utilities.contentsWithRegexpFromFolder(dosage_folder, dosage_pattern)
    files = sorted([os.path.join(dosage_folder, x) for x in files])
    for f in files:
        for line in model_training_file_geno_lines(f, snp_annotation, snps=snps):
            yield line

def model_training_geno_by_chromosome(dosage_folder, dosage_pattern, gtex_snp_file, gtex_release_version, snps=None):
    buffer = []
    last_chr = None

    snp_annotation = Helpers.gtex_snp(gtex_snp_file, gtex_release_version)

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
    for line in model_training_folder_geno_lines(dosage_folder, dosage_pattern, snp_annotation, snps):
        chromosome = line[GF.CHROMOSOME]

        if last_chr is None: last_chr = chromosome

        if last_chr != chromosome:
            yield _buffer_to_data(buffer)
            buffer = []

        last_chr = chromosome
        buffer.append(line)

    yield _buffer_to_data(buffer)


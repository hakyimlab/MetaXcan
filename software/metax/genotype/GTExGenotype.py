import gzip
import logging

import numpy

from Genotype import GF
from .. import Utilities
from .. import Exceptions
from ..misc import KeyedDataSource


def parse_gtex_variant(variant):
    comps = variant.split("_")
    return comps[0:4]

def gtex_geno_header(gtex_file):
    with gzip.open(gtex_file) as file:
        header = file.readline().strip().split()
    return header

def gtex_geno_lines(gtex_file, gtex_snp_file, snps=None, gtex_release_version=None):
    logging.log(9, "Loading GTEx snp file")
    #TODO: change to something more flexible to support V7 naming

    if not gtex_release_version:
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "VariantID", "RS_ID_dbSNP142_CHG37p13", numeric=False)
    elif gtex_release_version.lower() == "v7":
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "variant_id", "rs_id_dbSNP147_GRCh37p13", numeric=False)
    elif gtex_release_version.lower() == "v8":
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "variant_id", "rs_id_dbSNP150_GRCh38p7", numeric=False)
    else:
        raise Exceptions.InvalidArguments("Unsupported GTEx relkease version")

    logging.log(9, "Processing GTEx geno")
    with gzip.open(gtex_file) as file:
        for i,line in enumerate(file):
            if i==0:continue #skip header. This line is not needed but for conceptual ease of mind
            comps = line.strip().split()
            variant = comps[0]

            if not variant in gtex_snp:
                continue

            rsid = gtex_snp[variant]
            if snps and not rsid in snps:
                continue

            data = parse_gtex_variant(variant)
            dosage = numpy.array(comps[1:], dtype=numpy.float64)
            frequency = numpy.mean(dosage)/2

            yield [rsid] + data + [frequency] + list(dosage)

def gtex_geno_by_chromosome(gtex_file, gtex_snp_file, snps=None, gtex_release_version=None):
    buffer = []
    last_chr = None

    def _buffer_to_data(buffer):
        _data = zip(*buffer)

        _metadata = _data[0:GF.FIRST_DOSAGE]
        rsids = _metadata[GF.RSID]

        chr_ = _metadata[GF.CHROMOSOME][0].replace("chr", "")
        chromosome = int(chr_) # I know I am gonna regret this cast UPDATE: I did regret it!
        _metadata = zip(*_metadata)
        metadata = Utilities.to_dataframe(_metadata, ["rsid", "chromosome", "position", "ref_allele", "alt_allele", "frequency"], to_numeric="ignore")

        dosage = zip(*_data[GF.FIRST_DOSAGE:])
        dosage_data = {}
        #TODO: 142_snps are not unique, several rows go into the same value, improve this case handling
        for i in xrange(0, len(rsids)):
            rsid = rsids[i]
            if rsid in dosage_data:
                ind = metadata[metadata.rsid == rsid].index.tolist()[0]
                metadata = metadata.drop(ind)

            dosage_data[rsid] = dosage[i]
        metadata["number"] = range(0, len(metadata))
        metadata = metadata.set_index("number")

        return chromosome, metadata, dosage_data

    logging.log(8, "Starting to process lines")
    for line in gtex_geno_lines(gtex_file, gtex_snp_file, snps, gtex_release_version):

        chromosome = line[GF.CHROMOSOME]
        if last_chr is None: last_chr = chromosome

        if last_chr != chromosome:
            yield _buffer_to_data(buffer)
            buffer = []

        last_chr = chromosome
        buffer.append(line)

    yield _buffer_to_data(buffer)
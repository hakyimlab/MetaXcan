
from cyvcf2 import VCF
import logging
import pandas
import numpy

from metax.misc import Genomics

def vcf_file_geno_lines(path, mode="genotyped", variant_mapping=None, whitelist=None, skip_palindromic=False, liftover_conversion=None):
    logging.log(9, "Processing vcf %s", path)
    vcf_reader = VCF(path)

    is_dict_mapping = variant_mapping is not None and type(variant_mapping) == dict

    for variant in vcf_reader:
        chr = variant.CHROM
        pos = variant.POS
        variant_id = variant.ID
        ref = variant.REF
        alts = variant.ALT

        if liftover_conversion:
            chr_, pos_ = chr, pos
            chr, pos = liftover_conversion(chr, pos)
            if chr == "NA" or pos == "NA":
                continue

        if mode == "genotyped":
            for a,alt in enumerate(alts):
                if skip_palindromic and Genomics.is_palindromic(ref, alt):
                    continue

                _varid, variant_id = Genomics.maybe_map_variant(variant_id, chr, pos, ref, alt, variant_mapping, is_dict_mapping)
                if variant_id is None: continue

                if whitelist and variant_id not in whitelist:
                    continue

                d = []
                for sample in variant.genotypes:
                    d_ = (sample[0] == a+1) + (sample[1] == a+1)
                    d.append(d_)
                f = numpy.mean(numpy.array(d,dtype=numpy.int32))/2
                yield (variant_id, chr, pos, ref, alt, f) + tuple(d)

        elif mode == "imputed":
            if len(alts) > 1:
                logging.log("VCF imputed mode doesn't support multiple ALTs, skipping %s", variant_id)
                continue

            alt = alts[0]
            if skip_palindromic and Genomics.is_palindromic(ref, alt):
                continue

            _varid, variant_id = Genomics.maybe_map_variant(variant_id, chr, pos, ref, alt, variant_mapping, is_dict_mapping)
            if variant_id is None: continue

            if whitelist and variant_id not in whitelist:
                continue
            
            try:
                d = numpy.apply_along_axis(lambda x: x[0], 1, variant.format("DS"))
                f = numpy.mean(numpy.array(d)) / 2
                yield (variant_id, chr, pos, ref, alt, f) + tuple(d)
            except KeyError:
                yield RuntimeError("Missing DS field when vcf mode is imputed")
        else:
            yield RuntimeError(f"Unsupported vcf mode = {mode}")


def vcf_files_geno_lines(files, mode="genotyped", variant_mapping=None, whitelist=None, skip_palindromic=False, liftover_conversion=None):
    logging.log(9, "Processing vcfs")
    for file in files:
        for l in vcf_file_geno_lines(file, mode=mode, variant_mapping=variant_mapping,
                    whitelist=whitelist, skip_palindromic=skip_palindromic,
                    liftover_conversion=liftover_conversion):
            yield l

def get_samples(path):
    vcf_reader = VCF(path)
    ids = [(x, x) for x in vcf_reader.samples]
    ids = pandas.DataFrame(ids, columns=["FID", "IID"])
    return ids

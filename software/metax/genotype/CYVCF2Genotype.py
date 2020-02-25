
from cyvcf2 import VCF
import logging
import pandas
import numpy


def maybe_map_variant(varid, chr, pos, ref, alt, variant_mapping, is_dict_mapping):
    _varid = varid

    if variant_mapping:
        if is_dict_mapping:
            if not varid in variant_mapping:
                varid = None
            else:
                varid = variant_mapping[varid]
        else:
            varid = variant_mapping(chr, pos, ref, alt)
    return _varid, varid

def vcf_file_geno_lines(path, mode="genotyped", variant_mapping=None, whitelist=None):
    logging.log(9, "Processing vcf %s", path)
    vcf_reader = VCF(path)

    is_dict_mapping = variant_mapping is not None and type(variant_mapping) == dict

    for variant in vcf_reader:
        chr = variant.CHROM
        pos = variant.POS
        variant_id = variant.ID
        ref = variant.REF
        alts = variant.ALT

        if whitelist and variant_id not in whitelist:
            continue

        if mode == "genotyped":
            for a,alt in enumerate(alts):
                _varid, variant_id = maybe_map_variant(variant_id, chr, pos, ref, alt, variant_mapping, is_dict_mapping)
                if variant_id is None: continue

                d = []
                for sample in variant.genotypes:
                    d_ = (sample[0] == a+1) + (sample[1] == a+1)
                    d.append(d_)
                f = numpy.mean(numpy.array(d,dtype=numpy.int32))/2
                yield (variant_id, chr, pos, ref, alt, f) + tuple(d)

        elif mode == "imputed":
            if len(alts) > 1:
                raise RuntimeError("VCF imputed mode doesn't support multiple ALTs")
            alt = alts[0]
            _varid, variant_id = maybe_map_variant(variant_id, chr, pos, ref, alt, variant_mapping, is_dict_mapping)
            if variant_id is None: continue

            d = []
            for sample in variant.genotypes:
                d_ = sample[0] + sample[1]
                d.append(d_)
            f = numpy.mean(numpy.array(d, dtype=numpy.float64)) / 2
            yield (variant_id, chr, pos, ref, alt, f) + tuple(d)
        else:
            raise  RuntimeError("Unsupported vcf mode")


def vcf_files_geno_lines(files, mode="genotyped", variant_mapping=None, whitelist=None):
    logging.log(9, "Processing vcfs")
    for file in files:
        for l in vcf_file_geno_lines(file, mode, variant_mapping, whitelist):
            yield l

def get_samples(path):
    vcf_reader = VCF(path)
    ids = [(x, x) for x in vcf_reader.samples]
    ids = pandas.DataFrame(ids, columns=["FID", "IID"])
    return ids
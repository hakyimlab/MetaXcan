import bgen_reader
import numpy
import pandas
import logging

from ..misc import Genomics

def bgen_file_geno_lines(file, variant_mapping = None, force_colon = False, use_rsid=False, whitelist=None, skip_palindromic=False, liftover_conversion=None):
    logging.log(9, "Processing bgen %s", file)
    bgen = bgen_reader.read_bgen(file)
    variants = bgen["variants"]
    dict_mapping = variant_mapping is not None and type(variant_mapping) == dict
    for variant in variants.itertuples():
        if use_rsid:
            varid = variant.rsid
        else:
            varid = variant.id

        if force_colon:
            varid = varid.replace("_", ":")

        alleles = variant.allele_ids.split(",")
        if len(alleles) > 2:
            logging.info("variant %s is multiallelic, skipping", varid)
            continue
        allele_0, allele_1 = alleles[0], alleles[1]
        if skip_palindromic and Genomics.is_palindromic(allele_0, allele_1):
            continue

        pos = variant.pos
        chr = variant.chrom
        if liftover_conversion:
            chr_, pos_ = chr, pos
            chr, pos = liftover_conversion(chr, pos)
            if chr == "NA" or pos == "NA":
                continue

        if variant_mapping:
            if dict_mapping:
                if not varid in variant_mapping:
                    continue
                else:
                    varid_ = varid
                    varid = variant_mapping[varid]
            else:
                raise RuntimeError("BGEN code doessn't support variant mapping through a method")
        # subtlety: even though we replace the variant id,
        # the alleles in the genotype might be swapped respect the variant in the mapping
        # You should verify if you must match it


        if whitelist and not varid in whitelist:
            continue

        v = bgen["genotype"][variant.Index].compute()
        if v["phased"]:
            d = numpy.apply_along_axis(lambda x: x[1] + x[3], 1, numpy.array(v["probs"], dtype=numpy.float))
        else:
            d = numpy.apply_along_axis(lambda x: x[1] + x[2] * 2, 1, numpy.array(v["probs"], dtype=numpy.float))

        #e = bgen_reader.allele_expectation(bgen, variant.Index)
        #d2 = bgen_reader.compute_dosage(e, alt=1)

        yield (varid, chr, pos, allele_0, allele_1, numpy.mean(d)/2) + tuple(d)

def bgen_files_geno_lines(files, variant_mapping = None, force_colon = False, use_rsid=False, whitelist=None, skip_palindromic=False, liftover_conversion=None):
    logging.log(9, "Processing bgens")
    for file in files:
        for l in bgen_file_geno_lines(file, variant_mapping=variant_mapping, force_colon=force_colon, use_rsid=use_rsid, whitelist=whitelist, skip_palindromic=skip_palindromic, liftover_conversion=liftover_conversion):
            yield l

def get_samples(path):
    logging.info("Opening bgen to get samples")
    bgen = bgen_reader.read_bgen(path, verbose=False)
    samples = bgen["samples"].values
    samples = pandas.DataFrame({"FID":samples, "IID":samples})[["FID", "IID"]]
    return samples
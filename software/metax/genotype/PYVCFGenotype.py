import vcf
import logging
import pandas
import numpy

def int_treatment(x, allele_index):
    try:
        d_ = (int(x) == allele_index + 1)
    except:
        d_ = 0
    return d_

def float_treatment(x):
    try:
        d_ = float(x)
    except:
        d_ =0
    return d_

def vcf_file_geno_lines(path, mode="genotyped", whitelist=None):
    vcf_reader = vcf.Reader(filename=path)
    for record in vcf_reader:
        chr = record.CHROM
        pos = record.POS
        variant_id = record.ID
        ref = record.REF
        alts = record.ALT

        if whitelist and not variant_id in whitelist:
            continue

        if mode == "genotyped":
            for a,alt in enumerate(alts):
                d = []
                for sample in record.samples:
                    d_ = 0
                    d_ += int_treatment(sample.gt_alleles[0], a)
                    d_ += int_treatment(sample.gt_alleles[1], a)
                    d.append(d_)
                f = numpy.mean(numpy.array(d,dtype=numpy.int32))/2
                yield (variant_id, chr, pos, ref, alt, f) + tuple(d)
        elif mode == "imputed":
            if len(alts) > 1:
                raise RuntimeError("VCF imputed mode doesn't support multipl ALTs")
            d = []
            for sample in record.samples:
                d_ = 0
                d_ += float_treatment(sample.gt_alleles[0])
                d_ += float_treatment(sample.gt_alleles[1])
                d.append(d_)
            f = numpy.mean(numpy.array(d, dtype=numpy.float64)) / 2
            yield (variant_id, chr, pos, ref, alt, f) + tuple(d)
        else:
            raise  RuntimeError("Unsupported vcf mode")



def vcf_files_geno_lines(files, mode="genotyped", whitelist=None):
    logging.log(9, "Processing vcfs")
    for file in files:
        for l in vcf_file_geno_lines(file, mode, whitelist):
            yield l

def get_samples(path):
    vcf_reader = vcf.Reader(filename=path)
    r = next(vcf_reader)
    ids = [(x.sample, x.sample) for x in r.samples]
    ids = pandas.DataFrame(ids, columns=["FID", "IID"])
    return ids
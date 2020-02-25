

def map_on_the_fly(mapping, format, chromosome, position, ref_allele, alt_allele):
    v = format.format(chromosome, position, ref_allele, alt_allele)
    r = mapping[v] if v in mapping else None
    return r
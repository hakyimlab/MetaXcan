

def coordinate_format(format, chromosome, position, ref_allele, alt_allele):
    return format.format(chromosome, position, ref_allele, alt_allele)

def map_on_the_fly(mapping, format, chromosome, position, ref_allele, alt_allele):
    v = format.format(chromosome, position, ref_allele, alt_allele)
    r = mapping[v] if v in mapping else None
    return r

def is_palindromic(ref_allele, alt_allele):
    if ref_allele == "C" and alt_allele == "G":
        return True
    elif ref_allele == "G" and alt_allele == "C":
        return True
    elif ref_allele == "A" and alt_allele == "T":
        return True
    elif ref_allele == "T" and alt_allele == "A":
        return True
    return False
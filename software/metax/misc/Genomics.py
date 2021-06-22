import pyliftover
import logging

BASEPAIR = {
    "A": "T",
    "C": "G", 
    "G": "C",
    "T": "A"
}

def try_reverse_compelement(str_):
    '''
    If reverse complement fails, return empty string
    '''
    res = ''
    for i in str_:
        if i in BASEPAIR:
            res = BASEPAIR[i] + res
        else:
            return ''
    return res

def coordinate_format(checklist, format, chromosome, position, ref_allele, alt_allele):
    r = None
    # direct formatting
    v = format.format(chromosome, position, ref_allele, alt_allele)
    # flipped formatting
    v_ = format.format(chromosome, position, alt_allele, ref_allele)
    # complement formatting
    c = format.format(chromosome, position, try_reverse_compelement(ref_allele), try_reverse_compelement(alt_allele))
    # flipped complement formatting
    c_ = format.format(chromosome, position, try_reverse_compelement(alt_allele), try_reverse_compelement(ref_allele))

    if v in checklist:
        r = v
    elif v_ in checklist:
        r = v_
    elif c in checklist:
        r = c
    elif c_ in checklist:
        r = c_


    return r

def map_on_the_fly(mapping, format, chromosome, position, ref_allele, alt_allele):
    r = None
    v = format.format(chromosome, position, ref_allele, alt_allele)

    if v in mapping:
        r = mapping[v]

    elif len(ref_allele) == 1 and len(alt_allele) == 1:
        #try swapped
        v_ = format.format(chromosome, position, alt_allele, ref_allele)
        #try complement
        c = format.format(chromosome, position, try_reverse_compelement(ref_allele), try_reverse_compelement(alt_allele))
        # try swapped complement
        c_ = format.format(chromosome, position, try_reverse_compelement(alt_allele), try_reverse_compelement(ref_allele))

        if v_ in mapping:
            r = mapping[v_]

        elif c in mapping:              
            r = mapping[c]
        
        elif c_ in mapping:              
            r = mapping[c_]        

    return r

def is_palindromic(ref_allele, alt_allele):
    if try_reverse_compelement(ref_allele) == alt_allele:
        return True
    return False

def lift(liftover, chr, pos, zero_based_positions=False):
    _new_chromosome = "NA"
    _new_position = "NA"

    append_chromosome=False
    if not chr[0] == "c":
        append_chromosome =True
        chr = "chr"+chr
    try:
        pos = int(pos)
        if not zero_based_positions:
            pos = pos-1

        l_ = liftover.convert_coordinate(chr, pos)
        if l_:
            if len(l_) > 1:
                logging.warning("Liftover with more than one candidate: %s", t.variant_id)
            _new_chromosome = l_[0][0]
            _new_position = int(l_[0][1])

            if append_chromosome:
                _new_chromosome = _new_chromosome[3:]

            if not zero_based_positions:
                _new_position = _new_position+1
    except:
        pass
    return _new_chromosome, _new_position


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

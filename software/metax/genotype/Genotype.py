from ..misc import GWASAndModels

class GF(object):
    RSID=0
    CHROMOSOME=1
    POSITION=2
    REF_ALLELE=3
    ALT_ALLELE=4
    FREQUENCY=5
    FIRST_DOSAGE=6

def force_mapped_metadata(generator, sep):
    for line in generator:
        varid = line[GF.RSID]
        allele_0, allele_1 = line[GF.REF_ALLELE], line[GF.ALT_ALLELE]
        comps = varid.split(sep)
        f = line[GF.FREQUENCY]
        d = line[GF.FIRST_DOSAGE:]

        allele_swap, strand_swap = GWASAndModels.match_alleles(allele_0, allele_1, comps[2], comps[3])
        if not allele_swap or not strand_swap:
            continue

        pos = int(comps[1])
        chr = comps[0]
        if strand_swap == -1:
            allele_0 = comps[2]
            allele_1 = comps[3]
        if allele_swap == -1:
            allele_0, allele_1 = allele_1, allele_0
            f = 1-f
            d = tuple(map(lambda x:2-x, d))
        yield (varid, chr, pos, allele_0, allele_1, f) + d
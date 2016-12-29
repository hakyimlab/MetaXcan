def sample_dosage_data_1():
    #chromosome, rsid, position, non_effect_allele, effect_allele, frequency, allele_dosage
    s = [
        ("rs1", "chr1", 1, "T", "C", 0.1, [1, 0, 0, 0, 0]),
        ("rs2", "chr1", 5, "T", "C", 0.1, [0, 1, 0, 0, 0]),
        ("rs3", "chr1", 20, "G", "A", 0.2, [0, 0, 2, 0, 0]),
        ("rs4", "chr1", 30, "T", "C", 0.3, [1, 0, 1, 0, 1]),
        ("rs5", "chr1", 35, "A", "G", 0.7, [1, 0, 2, 1, 1]),  #
        ("rs6", "chr1", 42, "G", "A", 0.4, [1, 1, 1, 1, 0]),
        ("rs7", "chr1", 43, "C", "T", 0.7, [2, 0, 2, 2, 1]),
        ("rs8", "chr1", 50, "A", "G", 0.3, [0, 1, 0, 0, 2]),
        ("rs9", "chr1", 70, "G", "A", 0.3, [0, 0, 0, 1, 2]),
        ("rs10", "chr1", 55, "T", "G", 0.5, [1, 1, 1, 1, 1]),
        ("rs11", "chr1", 75, "C", "A", 1, [2, 2, 2, 2, 2]),
    ]

    return s

#similar as before but unordered
def sample_dosage_data_2():
    #rsid, chromosome,  position, non_effect_allele, effect_allele, frequency, allele_dosage
    s = [
        ("rsC", "chr1", 43, "C", "T", 0.7, [2, 0, 2, 2, 1]),
        ("rs1", "chr1", 1, "T", "C", 0.1, [1, 0, 0, 0, 0]),
        ("rs10", "chr1", 55, "T", "G", 0.5, [1, 1, 1, 1, 1]),
        ("rs4", "chr1", 30, "T", "C", 0.3, [1, 0, 1, 0, 1]),
        ("rs5", "chr1", 35, "A", "G", 0.7, [1, 0, 2, 1, 1]),  #
        ("rsB", "chr1", 42, "G", "A", 0.4, [1, 1, 1, 1, 0]),
        ("rs8", "chr1", 50, "A", "G", 0.3, [0, 1, 0, 0, 2]),
        ("rs9", "chr1", 70, "G", "A", 0.3, [0, 0, 0, 1, 2]),
        ("rs11", "chr1", 75, "C", "A", 1, [2, 2, 2, 2, 2]),
        ("rs2", "chr1", 5, "T", "C", 0.1, [0, 1, 0, 0, 0]),
        ("rsA", "chr1", 20, "G", "A", 0.2, [0, 0, 2, 0, 0]),
    ]

    return s

def sample_gwas_data_1():
    #SNP, CHROMOSOME, POSITION, NON_EFFECT_ALLELE, EFFECT_ALLELE, ZSCORE
    g = [
        ("rs1666", "chr1", 0, "G", "A", 0.3),
        ("rs1", "chr1", 1, "T", "C", -0.2),
        ("rs2", "chr1", 5, "T", "C", 0.5),
        ("rs3", "chr1", 20, "A", "G", 1.3),
        ("rs4", "chr1", 30, "G", "A", -0.3),
        ("rs6", "chr1", 42, "A", "G", 2.9),
        ("rs7", "chr1", 43, "C", "T", 4.35),
        ("rs7666", "chr1", 45, "G", "A", 1.3),
        ("rs8", "chr1", 50, "G", "A", 0.09),
        ("rs9", "chr1", 70, "G", "A", 0.09),
    ]
    return g

def sample_gwas_data_1_e():
    #SNP, CHROMOSOME, POSITION, NON_EFFECT_ALLELE, EFFECT_ALLELE, ZSCORE
    g = [
        ("rs1666", "chr1", 0, "G", "A", 0.3, 0, "a"),
        ("rs1", "chr1", 1, "T", "C", -0.2, 1, "b"),
        ("rs2", "chr1", 5, "T", "C", 0.5, 2, "c"),
        ("rs3", "chr1", 20, "A", "G", 1.3, 3, "d"),
        ("rs4", "chr1", 30,"G", "A", -0.3, 4, "e"),
        ("rs6", "chr1", 42, "A", "G", 2.9, 5, "f"),
        ("rs7", "chr1", 43, "C", "T", 4.35, 6, "g"),
        ("rs7666", "chr1", 45, "G", "A", 1.3, 7, "h"),
        ("rs8", "chr1", 50, "G", "A", 0.09, 8, "i"),
        ("rs9", "chr1", 70, "G", "A", 0.09, 9, "j"),
    ]
    return g

def sample_gwas_data_2():
    #CHROMOSOME, SNP, NON_EFFECT_ALLELE, EFFECT_ALLELE, ZSCORE
    g = [
        ("rsC", "chr1", None, "C", "T", 4.35),
        ("rs1666", "chr1", None, "G", "A", 0.3),
        ("rs1", "chr1", None, "T", "C", -0.2),
        ("rs2", "chr1", None, "T", "C", 1.3),
        ("rs4", "chr1", None, "G", "A", -0.3),
        ("rsB", "chr1", None, "A", "G", 2.9),
        ("rsA", "chr1", None, "A", "G", 1.3),
        ("rs7666", "chr1", None, "G", "A", 1.3),
        ("rs8", "chr1", None, "G", "A", 0.09),
        ("rs9", "chr1", None, "G", "A", 0.09),
    ]
    return g

# same as 1, but unordered
def sample_gwas_data_3():
    #CHROMOSOME, SNP, NON_EFFECT_ALLELE, EFFECT_ALLELE, ZSCORE
    g = [
        ("rs6", "chr1", None, "A", "G", 2.9),
        ("rs1666", "chr1", None, "G", "A", 0.3),
        ("rs2", "chr1", None, "T", "C", 1.3),
        ("rs4", "chr1", None, "G", "A", -0.3),
        ("rs7", "chr1", None, "C", "T", 4.35),
        ("rs3", "chr1", None, "A", "G", 1.3),
        ("rs7666", "chr1", None, "G", "A", 1.3),
        ("rs8", "chr1", None, "G", "A", 0.09),
        ("rs9", "chr1", None, "G", "A", 0.09),
        ("rs1", "chr1", None, "T", "C", -0.2),
    ]
    return g
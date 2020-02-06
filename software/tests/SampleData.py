import pandas

def _dataframe(names, data):
    columns = list(zip(*data))
    d = {names[i]:columns[i] for i in range(0,len(names))}
    d = pandas.DataFrame(d)
    d = d[names]
    return d

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

def sample_gwas_data_4():
    #SNP, CHROMOSOME, POSITION, NON_EFFECT_ALLELE, EFFECT_ALLELE, ZSCORE, BETA
    g = [
        ("rs1666", "chr1", 0, "G", "A", 0.3, 0.16),
        ("rs1", "chr1", 1, "T", "C", -0.2, -0.11),
        ("rs2", "chr1", 5, "T", "C", 0.5, 0.29),
        ("rs3", "chr1", 20, "A", "G", 1.3, 1),
        ("rs4", "chr1", 30, "G", "A", -0.3, -0.1),
        ("rs6", "chr1", 42, "A", "G", 2.9, 1),
        ("rs7", "chr1", 43, "C", "T", 4.35, 3),
        ("rs7666", "chr1", 45, "G", "A", 1.3, 5),
        ("rs8", "chr1", 50, "G", "A", 0.09, 1),
        ("rs9", "chr1", 70, "G", "A", 0.09, 1),
        ("rs100", "chr22", 1, "G", "A", 0.09, 0.01),
        ("rs101", "chr22", 2, "G", "A", 0.09, 0.01),
        ("rs102", "chr22", 3, "G", "A", 0.09, 0.01),
        ("rs202", "chr22", 5, "G", "A", 0.09, 0.01),
        ("rs401", "chr22", 10, "C", "T", 0.09, 0.01),
        ("rs402", "chr22", 11, "C", "T", 0.09, 0.01),
    ]
    return g

def dataframe_from_gwas(gwas):
    names = ["snp", "chromosome", "position", "non_effect_allele", "effect_allele", "zscore"]
    if len(gwas[0]) > 6: names.append("beta")
    return _dataframe(names, gwas)

def sample_weights_1():
    #rsid, gene, weight, ref_allele, eff_allele)
    w = [
        ("rs1666", "A", 0.2, "G", "A"),
        ("rs1", "A", 0.1, "T", "C"),
        ("rs2", "A", 0.4, "C", "T"),
        ("rs2666", "A", -0.5, "T", "C"),
        ("rs3", "B", 0.2, "A", "G"),
        ("rs6", "B", -0.2, "T", "C"),
        ("rs7", "B", 0.4, "C", "T"),
        ("rs7666", "B", -0.1, "G", "A"),
        ("rs8", "B", -0.1, "G", "A"),
        ("rs9", "B", 0.1, "G", "A"),
        ("rs100", "C", 0.5, "G", "A"),
        ("rs101", "C", 0.2, "G", "A"),
        ("rs102", "C", 0.3, "C", "T"),
        ("rs201", "D", 0.3, "C", "T"),
        ("rs202", "D", 0.4, "A", "G"),
        ("rs301", "E", 0.9, "C", "T"),
        ("rs401", "F", 0.2, "C", "T"),
        ("rs402", "F", 0.2, "C", "T"),
        ("rs402", "G", 0.2, "C", "T"),
    ]
    return w

def sample_extra_1():
    #gene, genename, n.snps.in.model, pred.perf.R2, pred.perf.pval, pred.perf.qval
    e = [
        ("A", "gene1", 4, 0.4, 0.04, 0.04),
        ("B", "gene2", 6, 0.6, 0.06, 0.06),
        ("C", "gene3", 3, 0.3, 0.03, 0.03),
        ("D", "gene4", 2, 0.2, 0.02, 0.02),
        ("E", "gene5", 1, 0.1, 0.01, 0.01),
        ("F", "gene6", 2, 0.22, 0.022, 0.022),
        ("G", "gene6", 1, 0.11, 0.011, 0.011),
    ]
    return e

def sample_weights_2():
    w = [
        ("rs1", "A", 0.2, "C", "T"),
        ("rs2", "A", 0.1, "A", "G"),
        ("rs3", "A", 0.05, "G", "A"),
        ("rs4", "B", 0.4, "T", "C"),
        ("rs5", "B", 0.3, "C", "T"),
        ("rs6", "C", 0.5, "T", "C"),
        ("rs1", "D", 0.6, "T", "C")
    ]
    return w

def sample_extra_2():
    e =[
        ("A", "gene1", 3, 0.9, 0.09, 0.091),
        ("B", "gene2", 2, 0.8, 0.08, 0.081),
        ("C", "gene3", 1, 0.7, 0.07, 0.071),
        ("D", "gene4", 1, 0.6, 0.06, 0.061),
    ]
    return e

def dataframe_from_weights(weights):
    names = ["rsid", "gene", "weight", "non_effect_allele", "effect_allele"]
    return _dataframe(names, weights)

def dataframe_from_extra(extra):
    names = ["gene", "gene_name", "n_snps_in_model", "pred_perf_r2", "pred_perf_pval", "pred_perf_qval"]
    return _dataframe(names, extra)

def sample_covariance_s_1():
    #gene, rsid1, rsid2, value
    c = [
        ("A", "rs1666", "rs1666", 0.25),
        ("A", "rs1666", "rs1", -0.0833333333333),
        ("A", "rs1666", "rs2", 0.166666666667),
        ("A", "rs1666", "rs2666", -0.166666666667),
        ("A", "rs1", "rs1", 0.25),
        ("A", "rs1", "rs2", 0.166666666667),
        ("A", "rs1", "rs2666", -0.166666666667),
        ("A", "rs2", "rs2", 0.333333333333),
        ("A", "rs2", "rs2666", -0.333333333333),
        ("A", "rs2666", "rs2666", 0.333333333333),
        ("B", "rs3", "rs3", 0.25),
        ("B", "rs3", "rs6", -0.0833333333333),
        ("B", "rs3", "rs7", 0.166666666667),
        ("B", "rs3", "rs7666", -0.166666666667),
        ("B", "rs3", "rs8", -0.166666666667),
        ("B", "rs3", "rs9", 0.166666666667),
        ("B", "rs6", "rs6", 0.25),
        ("B", "rs6", "rs7", 0.166666666667),
        ("B", "rs6", "rs7666", -0.166666666667),
        ("B", "rs6", "rs8", -0.166666666667),
        ("B", "rs6", "rs9", -0.166666666667),
        ("B", "rs7", "rs7", 0.333333333333),
        ("B", "rs7", "rs7666", -0.333333333333),
        ("B", "rs7", "rs8", -0.333333333333),
        ("B", "rs7", "rs9", 0.0),
        ("B", "rs7666", "rs7666", 0.333333333333),
        ("B", "rs7666", "rs8", 0.333333333333),
        ("B", "rs7666", "rs9", 0.0),
        ("B", "rs8", "rs8", 0.333333333333),
        ("B", "rs8", "rs9", 0.0),
        ("B", "rs9", "rs9", 0.333333333333),
        ("C", "rs100", "rs100", None),
        ("C", "rs100", "rs101", None),
        ("C", "rs100", "rs102", None),
        ("C", "rs101", "rs101", 0.333333333333),
        ("C", "rs101", "rs102", 0.0),
        ("C", "rs102", "rs102", 0.333333333333),
        ("F", "rs401", "rs401", 0.25),
        ("F", "rs401", "rs402", -0.25),
        ("F", "rs402", "rs402", 0.25),
        ("G", "rs402", "rs402", 0.0),
    ]
    return c

def dataframe_from_covariance(c):
    names = ["GENE", "RSID1", "RSID2", "VALUE"]
    return _dataframe(names, c)

def feature_set_1():
    #4 individuals, 2 genes
    c = [
        (0.0, 1.0),
        (0.1, 0.5),
        (0.3, 1.2),
        (0.0, 0.9)
    ]
    c = _dataframe(["a", "b"], c)
    return c

def feature_set_2():
    #4 individuals, 5 genes
    c = [
        (0.0, 1.0, 3.0, 0.5, 0.31),
        (0.1, 1.0, 0.0, 0.7, 0.32),
        (0.0, 1.2, 0.2, 0.6, 0.29),
        (0.1, 1.0, 0.3, 0.8, 0.3),
    ]
    c = _dataframe(["a", "b", "c", "d", "e"], c)
    return c

def feature_set_3():
    #4 individuals, 3 genes
    c = [
        (0.0, 1.0, 0.9),
        (0.1, 0.75, 0.9),
        (0.2, 1.2, 0.8),
        (0.1, 0.95, 1.1),
    ]
    c = _dataframe(["a", "b", "d"], c)
    return c

def set_of_feature_sets():
    sets = {}
    sets["1"] = feature_set_1()
    sets["2"] = feature_set_2()
    sets["3"] = feature_set_3()
    return sets

#!/usr/bin/env python
import gzip
import os
import string
import numpy
import pandas
from sqlalchemy import create_engine

#for reproducibility
numpy.random.seed(100)

sample_size = 100


def to_line(x):
    return "{}\n".format("\t".join([str(_x) for _x in x]))

def variant_name(chrom, pos, a, b):
    return "{}_{}_{}_{}_b37".format(chrom, pos, a, b)

complement = {"A": ["G", "C"], "T": ["G", "C"], "C":["T","A"], "G":["A","T"]}

def snp_data(chrom, pos):
    pos = chrom * 100 + pos
    a = numpy.random.choice(complement.keys())
    b = numpy.random.choice(complement[a])
    variant = variant_name(chrom, pos, a, b)
    return [variant]+list(numpy.random.uniform(high=2, low=0, size=sample_size))

def build_data():
    data = []
    for chrom in xrange(1,23):
        for pos in xrange(1,100):
            snp = snp_data(chrom, pos)
            data.append(snp)
    return data

def write_geno(data, file):
    with gzip.open(file, "w") as f:
        header = ["Id"] + ["GTEX-{}A".format(i) for i in xrange(0, sample_size)]
        f.write(to_line(header))

        for x in data:
            if numpy.random.uniform() > 0.95:
                continue
            line = to_line(x)
            f.write(line)


def write_snps(data, file):
    with gzip.open(file, "w") as f:
        header= ["Chr", "Pos", "VariantID", "Ref_b37", "Alt", "RS_ID_dbSNP135_original_VCF", "RS_ID_dbSNP142_CHG37p13", "Num_alt_per_site"]
        f.write(to_line(header))

        for i,x in enumerate(data):
            variant = x[0]
            comps = variant.split("_")
            rsid_135 = "rs"+str(i)

            u = numpy.random.uniform()
            rsid_142 = "rs"+str(i if u < 0.99 else i+1)

            line = [comps[0], comps[1], variant, comps[2], comps[3], rsid_135, rsid_142, 1]
            f.write(to_line(line))


gene_names = list(string.ascii_uppercase)
gene_chromosomes = [numpy.random.randint(1, 3) for i in xrange(0, 3)] + \
                   [numpy.random.randint(11, 13) for i in xrange(4, len(gene_names))]
def _simulate_model(snps, n_genes):
    def _simulate_gene_model(data, gene, chromosome):
        d = data[data.Chr == chromosome]
        n_snps = numpy.random.randint(1, 20)
        d = d.sample(n_snps)
        d = d.rename(columns={"RS_ID_dbSNP142_CHG37p13":"rsid", "Ref_b37":"ref_allele", "Alt":"eff_allele"})
        d = d[["rsid", "ref_allele", "eff_allele"]]
        d["gene"] = gene
        d["weight"] = numpy.random.uniform(-4,4,n_snps)
        return d

    weights = pandas.DataFrame()
    extra = pandas.DataFrame()
    for i in xrange(0, n_genes):
        gene = gene_names[i]
        chromosome = gene_chromosomes[i]
        w = _simulate_gene_model(snps, gene, chromosome)
        weights = pandas.concat([weights, w])

        pred = numpy.random.uniform()

        ed = {"gene": [gene], "genename": [gene.lower()], "n.snps.in.model": [len(w)], "pred.perf.R2":[pred], "pred.perf.pval":[pred], "pred.perf.qval":[pred]}
        e = pandas.DataFrame(data=ed)
        extra = pandas.concat([extra, e])
    weights = weights.sort_values(by="gene")
    return extra,weights


def write_models(snp_data_path, target_folder):
    def _get_snp_data(snp_data_path):
        snp_data = pandas.read_table(snp_data_path)
        d = []
        for t in snp_data.itertuples():
            chance = numpy.random.uniform()
            if chance > 0.80:
                allele_1 = t.Alt
                allele_2 = numpy.random.choice(complement[allele_1])
                #keep old variant name even if this fictitious snp will become inconsistent
                #variant = variant_name(t.Chr, t.Pos, allele_1, allele_2)
                t = t._replace(Ref_b37=allele_1, Alt=allele_2)
                d.append(t[1:])
            else:
                d.append(t[1:])
        d = zip(*d)
        header = snp_data.columns.values
        d = {header[i]:d[i] for i in xrange(0,len(header))}
        d = pandas.DataFrame(data=d)
        return d
    snps = _get_snp_data(snp_data_path)

    names = ["model_sim_1.db", "model_sim_2.db"]
    for i,name in enumerate(names):
        extra, weights = _simulate_model(snps, 10+i)
        path = os.path.join(target_folder,name)
        engine = create_engine("sqlite:///"+path)

        with engine.connect() as conn, conn.begin():
            extra.to_sql("extra", conn, index=False,if_exists="replace")
            weights.to_sql("weights",conn, index=False,if_exists="replace")


if __name__ == "__main__":
    gtex_file = "_td/genotype/gtex_like.txt.gz"
    gtex_snp_file = "_td/genotype/gtex_snp.txt.gz"
    weight_folder = "_td/dbs_2"
    if not os.path.exists(weight_folder): os.makedirs(weight_folder)
    folder = os.path.split(gtex_file)[0]
    if not os.path.exists(folder): os.makedirs(folder)
    data = build_data()
    write_geno(data, gtex_file)
    write_snps(data, gtex_snp_file)
    write_models("_td/genotype/gtex_snp.txt.gz", weight_folder)






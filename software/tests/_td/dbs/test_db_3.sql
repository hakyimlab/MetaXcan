CREATE TABLE extra (gene TEXT, genename TEXT, `n.snps.in.model` INTEGER, `pred.perf.R2` DOUBLE, `pred.perf.pval` DOUBLE, `pred.perf.qval` DOUBLE);
CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER);
CREATE INDEX extra_gene ON extra (gene);
CREATE INDEX weights_gene ON weights (gene);
CREATE INDEX weights_rsid ON weights (rsid);
CREATE INDEX weights_rsid_gene ON weights (rsid, gene);

insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("A", "gene1", 5, "NA", "NA", "NA");
insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("B", "gene2", 2, "NA", "NA", "NA");

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs940550", "A", 0.2, "G", "C");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs6650104", "A", -0.1, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs6594028", "A", 0.15, "T", "C");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs659666", "A", 0.15, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs12082473", "A", 0.15, "C", "T");

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs3131971", "B", 0.5, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs61770173", "B", -0.3, "C", "A");
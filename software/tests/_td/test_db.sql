CREATE TABLE extra (gene TEXT, genename TEXT, R2 DOUBLE,  `n.snps` INTEGER, pval DOUBLE);
CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER);
CREATE INDEX extra_gene ON extra (gene);
CREATE INDEX weights_gene ON weights (gene);
CREATE INDEX weights_rsid ON weights (rsid);
CREATE INDEX weights_rsid_gene ON weights (rsid, gene);

insert into extra(gene, genename, R2, `n.snps`, pval) values ("A", "gene1", 0.9, 3, 0.09);
insert into extra(gene, genename, R2, `n.snps`, pval) values ("B", "gene2", 0.8, 2, 0.08);
insert into extra(gene, genename, R2, `n.snps`, pval) values ("C", "gene3", 0.7, 1, 0.07);
insert into extra(gene, genename, R2, `n.snps`, pval) values ("D", "gene4", 0.6, 1, 0.06);

insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs1", "A", 0.2, "C", "T", 0.1, 3, 1);
insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs2", "A", 0.1, "A", "G", 0.2, 3, 2);
insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs3", "A", 0.05, "G", "A", 0.3, 3, 3);

insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs4", "B", 0.4, "T", "C", 0.4, 2, 4);
insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs5", "B", 0.3, "C", "T", 0.5, 2, 5);

insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs6", "C", 0.5, "T", "C", 0.6, 1, 6);

insert into weights(rsid, gene, weight, ref_allele, eff_allele, pval, N, cis) values ("rs1", "D", 0.6, "T", "C", 0.7, 1, 7);
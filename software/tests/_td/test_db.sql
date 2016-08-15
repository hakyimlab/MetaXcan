CREATE TABLE extra (gene TEXT, genename TEXT, `n.snps.in.model` INTEGER, `pred.perf.R2` DOUBLE, `pred.perf.pval` DOUBLE, `pred.perf.qval` DOUBLE);
CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER);
CREATE INDEX extra_gene ON extra (gene);
CREATE INDEX weights_gene ON weights (gene);
CREATE INDEX weights_rsid ON weights (rsid);
CREATE INDEX weights_rsid_gene ON weights (rsid, gene);

insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("A", "gene1", 3, 0.9, 0.09, 0.091);
insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("B", "gene2", 2, 0.8, 0.08, 0.081);
insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("C", "gene3", 1, 0.7, 0.07, 0.071);
insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("D", "gene4", 1, 0.6, 0.06, 0.061);

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs1", "A", 0.2, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs2", "A", 0.1, "A", "G");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs3", "A", 0.05, "G", "A");

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs4", "B", 0.4, "T", "C");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs5", "B", 0.3, "C", "T");

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs6", "C", 0.5, "T", "C");

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs1", "D", 0.6, "T", "C");
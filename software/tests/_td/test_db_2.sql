--mimicking Muscle Skeletal
CREATE TABLE extra (gene TEXT, genename TEXT, `n.snps.in.model` INTEGER, `pred.perf.R2` DOUBLE, `pred.perf.pval` DOUBLE, `pred.perf.qval` DOUBLE);
CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER);
CREATE INDEX extra_gene ON extra (gene);
CREATE INDEX weights_gene ON weights (gene);
CREATE INDEX weights_rsid ON weights (rsid);
CREATE INDEX weights_rsid_gene ON weights (rsid, gene);

insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("ENSG00000153814.7", "JAZF1", 17, 0.0678427838866996, 5.2052946696116e-07, 6.14282057719623e-07);
insert into extra(gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`, `pred.perf.qval`) values ("ENSG00000178363.3", "CALML3", 3,  0.0120792285976449,	0.0368619107435424,	0.0177311824057577);


insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs245915" "ENSG00000153814.7" 0.0009325161, "G", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs245913" "ENSG00000153814.7" 0.0150080934, "A", "G");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs245909" "ENSG00000153814.7" 0.0323974467, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs245906" "ENSG00000153814.7" 0.0032068109, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs10486599" "ENSG00000153814.7" 0.0965340659, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs144012121" "ENSG00000153814.7" 0.2255904, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs117887801" "ENSG00000153814.7" 0.0182896251, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs542000" "ENSG00000153814.7" 0.0118210293, "A", "G");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs544632" "ENSG00000153814.7" 0.0452374923, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs498475" "ENSG00000153814.7" 0.0238776449, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs849327" "ENSG00000153814.7" 0.0137637664, "A", "G");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs849336" "ENSG00000153814.7" 0.0130848652, "A", "G");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs849335" "ENSG00000153814.7" 0.0133741745, "T", "C");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs1513272" "ENSG00000153814.7" 0.0200467059, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs849135" "ENSG00000153814.7" 0.0157369455, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs849134" "ENSG00000153814.7" 0.020611184, "A", "G");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs860262" "ENSG00000153814.7" 0.0222322401, "C", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs849133" "ENSG00000153814.7" 0.0232020046, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs1635852" "ENSG00000153814.7" 0.0825764319, "T", "C");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs864745" "ENSG00000153814.7" 0.0058809985, "T", "C");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs112751321" "ENSG00000153814.7" 0.101057192, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs144273091" "ENSG00000153814.7" -0.0570374216, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs117462481" "ENSG00000153814.7" 0.1467541295, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs149305679" "ENSG00000153814.7" 0.0203423858, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs643036" "ENSG00000153814.7" 0.0082011566, "C", "T");

insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs1937888", "ENSG00000178363.3" , 0.010741086, "C", "T");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs17155745", "ENSG00000178363.3", -0.007338589, "G", "A");
insert into weights(rsid, gene, weight, ref_allele, eff_allele) values ("rs62626328", "ENSG00000178363.3", 0.01166571, "A", "G");

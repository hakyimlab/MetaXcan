#!/usr/bin/env Rscript

df <- read.csv("./results/allele_stats.csv")
df$wdb_ref_allele <- as.character(df$wdb_ref_allele)
df$wdb_eff_allele <- as.character(df$wdb_eff_allele)
df$legend_ref_allele <- as.character(df$legend_ref_allele)
df$legend_eff_allele <- as.character(df$legend_eff_allele)
df$gwas_ref_allele <- as.character(df$gwas_ref_allele)
df$gwas_eff_allele <- as.character(df$gwas_eff_allele)

print("1")

df$flip <-  ifelse(df$wdb_ref_allele==df$legend_eff_allele & df$wdb_eff_allele==df$legend_ref_allele, 1, 0)
total_flip <- sum(df$flip, na.rm=TRUE)
cat("flip:", total_flip, "\n")

df$same <-  ifelse(df$wdb_ref_allele==df$legend_ref_allele & df$wdb_eff_allele==df$legend_eff_allele & !is.na(df$legend_ref_allele), 1, 0)
total_same <- sum(df$same, na.rm=TRUE)
cat("same:",total_same,"\n")

df$diff <-  ifelse((df$wdb_ref_allele!=df$legend_ref_allele | df$wdb_eff_allele!=df$legend_eff_allele) & !is.na(df$legend_ref_allele), 1, 0)
total_diff <- sum(df$diff, na.rm=TRUE)
cat("dif:",total_diff,"\n")

df$na <- ifelse( is.na(df$legend_ref_allele), 1, 0)
total_na = sum(df$na, na.rm=TRUE)
cat("na:",total_na,"\n")

print("2")

df$flip2 <-  ifelse(df$wdb_ref_allele==df$gwas_eff_allele & df$wdb_eff_allele==df$gwas_ref_allele, 1, 0)
total_flip2 <- sum(df$flip2, na.rm=TRUE)
cat("flip2:", total_flip2, "\n")

df$same2 <-  ifelse(df$wdb_ref_allele==df$gwas_ref_allele & df$wdb_eff_allele==df$gwas_eff_allele & !is.na(df$gwas_ref_allele), 1, 0)
total_same2 <- sum(df$same2, na.rm=TRUE)
cat("same2:",total_same2,"\n")

df$diff2 <-  ifelse((df$wdb_ref_allele!=df$gwas_ref_allele | df$wdb_eff_allele!=df$gwas_eff_allele) & !is.na(df$gwas_ref_allele), 1, 0)
total_diff2 <- sum(df$diff2, na.rm=TRUE)
cat("dif2:",total_diff2,"\n")

df$na2 <- ifelse( is.na(df$gwas_ref_allele), 1, 0)
total_na2 = sum(df$na2)
cat("na2:",total_na2,"\n")

print( length(df$wdb_eff_allele) )
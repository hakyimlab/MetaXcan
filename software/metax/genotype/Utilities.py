import GTExGenotype

def genotype_by_chromosome_from_args(args, snps=None):
    return GTExGenotype.gtex_geno_by_chromosome(args.gtex_genotype_file, args.gtex_snp_file, snps)
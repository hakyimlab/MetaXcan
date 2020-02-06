from . import GTExGenotype
from . import DosageGenotype
from . import ModelTrainingGenotype

from .. import Exceptions

def genotype_by_chromosome_from_args(args, snps=None):
    if args.gtex_genotype_file:
        return GTExGenotype.gtex_geno_by_chromosome(args.gtex_genotype_file, args.gtex_snp_file, snps, args.gtex_release_version, args.impute_to_mean)
    elif args.dosage_genotype_folder:
        return DosageGenotype.dosage_geno_by_chromosome(args.dosage_genotype_folder, args.dosage_genotype_pattern, snps)
    elif args.model_training_genotype_folder:
        return ModelTrainingGenotype.model_training_geno_by_chromosome(args.model_training_genotype_folder, args.model_training_genotype_pattern, args.gtex_snp_file, args.gtex_release_version, snps)
    else:
        raise Exceptions.InvalidArguments("Invalid genotype")
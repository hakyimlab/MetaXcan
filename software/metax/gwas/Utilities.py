import logging
import pandas
import numpy

import GWAS
from .. import Constants
from .. import Utilities as BUtilities

def add_gwas_arguments_to_parser(parser):
    parser.add_argument("--snp_column", help="Name of snp column", default=None)
    parser.add_argument("--effect_allele_column", help="Name of effect allele column", default=None)
    parser.add_argument("--non_effect_allele_column", help="Name of alternative (non effect) allele column", default=None)
    parser.add_argument("--chromosome_column", help="Name of chromosome column", default=None)
    parser.add_argument("--position_column", help="Name of base position column", default=None)
    parser.add_argument("--freq_column", help="Name of frequency column", default=None)
    parser.add_argument("--beta_column", help="Name of snp beta column", default=None)
    parser.add_argument("--se_column", help="Name of snp beta standard error column", default=None)
    parser.add_argument("--or_column", help="Name of snp Odds ratio column", default=None)
    parser.add_argument("--zscore_column", help="Name of snp Z-Score ratio column", default=None)
    parser.add_argument("--pvalue_column", help="Name of snp p-value column", default=None)

def add_gwas_format_json_to_parser(parser):
    parser.add_argument("--input_gwas_format_json",
                        help="File containing a json description of the gwas.")

def override_gwas_format_dict_from_parameters(dict, parameters):
    if parameters.snp_column: dict[GWAS.COLUMN_SNP] = parameters.snp_column

    if parameters.effect_allele_column: dict[GWAS.COLUMN_EFFECT_ALLELE] = parameters.eff_allele_column

    if parameters.non_effect_allele_column: dict[GWAS.COLUMN_NON_EFFECT_ALLELE] = parameters.non_eff_allele_column

    if parameters.chromosome_column: dict[GWAS.COLUMN_CHROMOSOME] = parameters.chromosome_column

    if parameters.position_column: dict[GWAS.COLUMN_POSITION] = parameters.position_column

    if parameters.beta_column: dict[GWAS.COLUMN_BETA] = parameters.beta_column

    if parameters.se_column: dict[GWAS.COLUMN_SE] = parameters.se_column

    if parameters.or_column: dict[GWAS.COLUMN_OR] = parameters.or_column

    if parameters.zscore_column: dict[GWAS.COLUMN_ZSCORE] = parameters.zscore_column

    if parameters.freq_column: dict[GWAS.COLUMN_FREQ] = parameters.freq_column

    if parameters.pvalue_column: dict[GWAS.COLUMN_PVALUE] = parameters.pvalue_column

def gwas_parameters_from_args(args):
    gwas_format = {}
    if args.input_gwas_format_json:
        logging.info("Reading GWAS input from json file: %s", args.input_gwas_format_json)
        gwas_format = BUtilities.load_json(args.input_gwas_format_json)

    logging.info("Processing GWAS command line parameters")
    override_gwas_format_dict_from_parameters(gwas_format, args)
    return gwas_format

def gwas_from_data(data, extra_columns=None):
    """"Data should be a list of tuples as in [(snp, chromosome, non_effect_allele, effect_allele, zscore)]"""

    if len(data):
        d = zip(*data)
        F = GWAS.GWASF
        rsid, chromosome, position, non_effect_allele, effect_allele, zscore = d[F.SNP], d[F.CHROMOSOME], d[F.POSITION], d[F.NON_EFFECT_ALLELE], d[F.EFFECT_ALLELE], d[F.ZSCORE]
    else:
        rsid, chromosome, position, non_effect_allele, effect_allele, zscore = [], [], [], [], [], []

    g = pandas.DataFrame({Constants.SNP:numpy.array(rsid, dtype=numpy.str),
                        Constants.CHROMOSOME:numpy.array(chromosome, dtype=numpy.str),
                        Constants.POSITION:numpy.array(position),
                        Constants.EFFECT_ALLELE:numpy.array(effect_allele, dtype=numpy.str),
                        Constants.NON_EFFECT_ALLELE:numpy.array(non_effect_allele, dtype=numpy.str),
                        Constants.ZSCORE:numpy.array(zscore)})
    if len(data) and extra_columns:
        for k,i in extra_columns:
            g[k] = numpy.array(d[i])
    return g
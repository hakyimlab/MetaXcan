import logging
import gzip
import pandas
import numpy

import GWAS
from .. import Exceptions
from .. import Constants
from .. import Utilities as BUtilities

def add_gwas_arguments_to_parser(parser):
    parser.add_argument("--snp_column", help="Name of -snp column- in GWAS input file", default="SNP")
    parser.add_argument("--effect_allele_column", help="Name of -effect allele column- in GWAS input file", default="A1")
    parser.add_argument("--non_effect_allele_column", help="Name of -non effect allele column- in GWAS input file", default="A2")
    parser.add_argument("--chromosome_column", help="Name of -chromosome column- in GWAS input file", default=None)
    parser.add_argument("--position_column", help="Name of -base position column- in GWAS input file", default=None)
    parser.add_argument("--freq_column", help="Name of -frequency column- in GWAS input file", default=None)
    parser.add_argument("--beta_column", help="Name of snp association's -beta column- in GWAS input file", default=None)
    parser.add_argument("--beta_sign_column", help="Name of snp association's -sign of beta column- in GWAS input file", default=None)
    parser.add_argument("--or_column", help="Name of snp association's -odds ratio column- in GWAS input file", default=None)
    parser.add_argument("--se_column", help="Name of snp association's -beta standard error- column in GWAS input file", default=None)
    parser.add_argument("--zscore_column", help="Name of snp association's -Z-Score ratio column- in GWAS input file", default=None)
    parser.add_argument("--pvalue_column", help="Name of snp association's -p-value column- in GWAS input file", default=None)

def add_gwas_format_json_to_parser(parser):
    parser.add_argument("--input_gwas_format_json",
                        help="File containing a json description of the gwas.")

def override_gwas_format_dict_from_parameters(dict, parameters):
    if parameters.snp_column: dict[GWAS.COLUMN_SNP] = parameters.snp_column

    if parameters.effect_allele_column: dict[GWAS.COLUMN_EFFECT_ALLELE] = parameters.effect_allele_column

    if parameters.non_effect_allele_column: dict[GWAS.COLUMN_NON_EFFECT_ALLELE] = parameters.non_effect_allele_column

    if hasattr(parameters, 'chromosome_column') and parameters.chromosome_column: dict[GWAS.COLUMN_CHROMOSOME] = parameters.chromosome_column

    if hasattr(parameters, 'position_column') and parameters.position_column: dict[GWAS.COLUMN_POSITION] = parameters.position_column

    if hasattr(parameters, 'freq_column') and parameters.freq_column: dict[GWAS.COLUMN_FREQ] = parameters.freq_column

    if parameters.beta_column: dict[GWAS.COLUMN_BETA] = parameters.beta_column

    if parameters.beta_sign_column: dict[GWAS.COLUMN_BETA_SIGN] = parameters.beta_sign_column

    if parameters.se_column: dict[GWAS.COLUMN_SE] = parameters.se_column

    if parameters.or_column: dict[GWAS.COLUMN_OR] = parameters.or_column

    if parameters.zscore_column: dict[GWAS.COLUMN_ZSCORE] = parameters.zscore_column

    if parameters.pvalue_column: dict[GWAS.COLUMN_PVALUE] = parameters.pvalue_column

def gwas_format_from_args(args):
    gwas_format = {}
    if hasattr(args, "input_gwas_format_json") and args.input_gwas_format_json:
        logging.info("Reading GWAS input from json file: %s", args.input_gwas_format_json)
        gwas_format = BUtilities.load_json(args.input_gwas_format_json)

    logging.info("Processing GWAS command line parameters")
    override_gwas_format_dict_from_parameters(gwas_format, args)
    return gwas_format

def gwas_filtered_source(path, snps=None, snp_column_name=None, skip_until_header=None, separator=None):
    s = {}
    o = gzip.open if ".gz" in path else open
    with o(path) as file:
        header = None
        if skip_until_header:
            for line in file:
                if skip_until_header in line:
                    header = skip_until_header
                    c = line.split(skip_until_header)
                    if len(c) > 1: header += c[1]
                    break
            if header is None: raise Exceptions.ReportableException("Did not find specified header")
        else:
            header = file.readline()

        header_comps = header.strip().split(separator)
        s = {c:[] for c in header_comps}
        index = -1
        if snp_column_name:
            if not snp_column_name in header_comps: raise Exceptions.ReportableException("Did not find snp colum name")
            index = header_comps.index(snp_column_name)

        header_count = {k:header_comps.count(k) for k in header_comps}
        if len(header_count) < len(header_comps):
            duplicated = [k for k,v in header_count.iteritems() if v>1]
            logging.log("The input GWAS has duplicated columns: %s, will only use the first one in each case", str(duplicated))

        for line in file:
            comps = line.strip().split(separator)
            if snps and not comps[index] in snps:
                continue

            # Load only the first column if in presence of duplicated columns. Yuck!
            sentinel=set()
            for i,c in enumerate(comps):
                comp = header_comps[i]
                if comp in sentinel: continue
                sentinel.add(comp)
                c = sanitize_component(c)
                s[comp].append(c)

        for c in header_comps:
            s[c] = numpy.array(pandas.to_numeric(s[c], errors='ignore'))

    return s

import re
non_en_number = re.compile("^[-\+]?[0-9]*,{1}[0-9]+([eE]{1}[-\+]?[0-9]+)?$")
def sanitize_component(c):
    if non_en_number.match(c): c = c.replace(",",".")
    if c == "NA": c = None
    if c == ".": c = None
    return c


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

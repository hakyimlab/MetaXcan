import logging
import re
import os
import pandas
import numpy

from . import GWAS
from .. import Exceptions
from .. import Constants
from .. import Utilities as BUtilities

def add_gwas_arguments_to_parser(parser):
    parser.add_argument("--snp_column", help="Name of -snp column- in GWAS input file", default="SNP")
    parser.add_argument("--effect_allele_column", help="Name of -effect allele column- in GWAS input file", default="A1")
    parser.add_argument("--non_effect_allele_column", help="Name of -non effect allele column- in GWAS input file", default="A2")
    parser.add_argument("--chromosome_column", help="Name of -chromosome column- in GWAS input file")
    parser.add_argument("--position_column", help="Name of -base position column- in GWAS input file")
    parser.add_argument("--freq_column", help="Name of -frequency column- in GWAS input file")
    parser.add_argument("--beta_column", help="Name of snp association's -beta column- in GWAS input file")
    parser.add_argument("--beta_sign_column", help="Name of snp association's -sign of beta column- in GWAS input file")
    parser.add_argument("--or_column", help="Name of snp association's -odds ratio column- in GWAS input file")
    parser.add_argument("--se_column", help="Name of snp association's -beta standard error- column in GWAS input file")
    parser.add_argument("--zscore_column", help="Name of snp association's -Z-Score ratio column- in GWAS input file")
    parser.add_argument("--pvalue_column", help="Name of snp association's -p-value column- in GWAS input file")
    parser.add_argument("--separator", help="Character or string separating fields in input file. Defaults to any whitespace.")

    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                        " Specify this option (string value) to identify a header up to which file contents should be skipped.")

    parser.add_argument("--handle_empty_columns", default=False, action="store_true",
                    help="Some files have empty columns, with values not even coded to missing. This instructs the parser to handle those lines."
                    "Be sur eyou want to use this.")

    parser.add_argument("--input_pvalue_fix", help="If input GWAS pvalues are too small to handle, replace with these significance level. Use -0- to disable this behaviour and discard problematic snps.", type=int, default=1e-50)

    parser.add_argument("--keep_non_rsid", help="Keep non-rsid snps", action="store_true")

    parser.add_argument("--snp_map_file", help="table specifying conversion between a particular set of snps and those in the models' reference")

    parser.add_argument("--split_column", help="Present for future compatibility.", nargs="+")

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

def gwas_from_data(data, extra_columns=None):
    """"Data should be a list of tuples as in [(snp, chromosome, non_effect_allele, effect_allele, zscore)]"""

    if len(data):
        d = list(zip(*data))
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

def load_plain_gwas_from_args(args):
    regexp = re.compile(args.gwas_file_pattern) if args.gwas_file_pattern else  None
    gwas_format = gwas_format_from_args(args)
    GWAS.validate_format_basic(gwas_format)
    GWAS.validate_format_for_strict(gwas_format)

    _l = lambda x: GWAS.load_gwas(x, gwas_format, skip_until_header=args.skip_until_header,
            separator=args.separator, handle_empty_columns=args.handle_empty_columns, input_pvalue_fix=args.input_pvalue_fix,
            keep_non_rsid=args.keep_non_rsid)
    if args.gwas_folder:
        names = BUtilities.contentsWithRegexpFromFolder(args.gwas_folder, regexp)
        names.sort()  # cosmetic, because different filesystems/OS yield folders in different order
        files = [os.path.join(args.gwas_folder,x) for x in names]
        files = [_l(x) for x in files]
        gwas = pandas.concat(files)
    elif args.gwas_file:
        gwas = _l(args.gwas_file)
    else:
        raise Exceptions.InvalidArguments("Provide either --gwas_file or (--gwas_folder with optional --gwas_file_pattern)")
    return gwas

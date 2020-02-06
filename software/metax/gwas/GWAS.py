import pandas #likely to go away for a streaming approach
import numpy
import logging
import gzip
import scipy.stats as stats

from . import GWASSpecialHandling

from .. import  Exceptions

from ..Constants import SNP
from ..Constants import EFFECT_ALLELE
from ..Constants import NON_EFFECT_ALLELE
from ..Constants import ZSCORE
from ..Constants import CHROMOSOME
from ..Constants import POSITION
from ..Constants import OR
from ..Constants import BETA
from ..Constants import BETA_SIGN
from ..Constants import SE
from ..Constants import PVALUE

COLUMN_SNP="column_snp"
COLUMN_EFFECT_ALLELE="column_effect_allele"
COLUMN_NON_EFFECT_ALLELE="column_non_effect_allele"
COLUMN_CHROMOSOME="column_chromosome"
COLUMN_POSITION="column_position"
COLUMN_FREQ="column_freq"
COLUMN_BETA="column_beta"
COLUMN_BETA_SIGN="column_beta_sign"
COLUMN_PVALUE="column_pvalue"
COLUMN_SE="column_se"
COLUMN_OR="column_or"
COLUMN_ZSCORE="column_zscore"

########################################################################################################################
# Format Management

#format of raw results
class GWASF(object):
    SNP=0
    CHROMOSOME=1
    POSITION=2
    NON_EFFECT_ALLELE=3
    EFFECT_ALLELE=4
    ZSCORE=5

def _f_snp(format): return format[COLUMN_SNP] if COLUMN_SNP in format else None
def _f_effect_allele_column(format): return format[COLUMN_EFFECT_ALLELE] if COLUMN_EFFECT_ALLELE in format else None
def _f_non_effect_allele_column(format): return format[COLUMN_NON_EFFECT_ALLELE] if COLUMN_NON_EFFECT_ALLELE in format else None
def _f_zscore(format): return format[COLUMN_ZSCORE] if COLUMN_ZSCORE in format else None
def _f_pvalue(format): return format[COLUMN_PVALUE] if COLUMN_PVALUE in format else None
def _f_beta(format): return format[COLUMN_BETA] if COLUMN_BETA in format else None
def _f_beta_sign(format): return format[COLUMN_BETA_SIGN] if COLUMN_BETA_SIGN in format else None
def _f_or(format): return format[COLUMN_OR] if COLUMN_OR in format else None
def _f_se(format): return format[COLUMN_SE] if COLUMN_SE in format else None

def validate_format_basic(format):
    if not _f_snp(format): raise Exceptions.InvalidArguments("Need to provide a SNP column")
    if not _f_effect_allele_column(format): raise Exceptions.InvalidArguments("Need to provide an -effect allele- column")
    if not _f_non_effect_allele_column(format): raise Exceptions.InvalidArguments("Need to provide a -non effect allele- column")

def validate_format_for_strict(format):
    ok = False
    if _f_zscore(format):
        ok = True
    elif _f_pvalue(format):
        if _f_beta(format) or _f_or(format) or _f_beta_sign(format):
            ok = True
        else:
            raise Exceptions.InvalidArguments("If providing pvalue, you must provide either -beta-, -sign of beta-, or -odd ratio-")
    elif _f_se(format):
        if _f_beta(format) or _f_or(format):
            ok = True
        else:
            raise Exceptions.InvalidArguments("If providing standard error, you must provide either -beta- or -odd ratio-")

    if not ok:
        raise Exceptions.InvalidArguments("Arguments missing. Either -zscore-, -pvalue and one in [beta,beta sign,or]-, or -standard error and one in [beta, or]- must be provided")


########################################################################################################################
# Load a gwas
def load_gwas(source, gwas_format, strict=True, separator=None, skip_until_header=False, snps=None, force_special_handling=False, handle_empty_columns=False, input_pvalue_fix=None, keep_non_rsid=False):
    """
    Attempts to read a GWAS summary statistics file, and load it into a uniform format,
    in a pandas dataframe.

    :param source: Either a string with path to file containing GWAS summary statistics, or a generator.
    :param gwas_format: dictionary specifying GWAS format column mapping.
        For example
    :return:
    """
    if force_special_handling or skip_until_header or snps:
        logging.info("Reading input gwas with special handling: %s", source)
        snp_column_name = gwas_format[COLUMN_SNP]
        d = GWASSpecialHandling.gwas_data_source(source, snps, snp_column_name, skip_until_header, separator, handle_empty_columns)
        d = pandas.DataFrame(d)
    else:
        logging.info("Reading input gwas: %s", source)
        if separator is None or separator == "ANY_WHITESPACE":
            separator = '\s+'
        d = pandas.read_table(source, separator)

    logging.info("Processing input gwas")
    d = _rename_columns(d, gwas_format)

    if not SNP in d:
        raise Exceptions.ReportableException("A valid SNP column name must be provided in the format")

    #keep only rsids
    if d.shape[0] > 0:
        d = d[~ d[SNP].isnull()]
        if not keep_non_rsid:
            d = d[d[SNP].str.contains("rs")]

    if strict:
        d = _enforce_numeric_columns(d)
        d = _ensure_columns(d, input_pvalue_fix)
        d = _keep_gwas_columns(d)
        if d.shape[0] >0 and numpy.any(~ numpy.isfinite(d[ZSCORE])):
            logging.warning("Some GWAS snp zscores are not finite.")

    return d

def _keep_gwas_columns(d):
    keep_columns = [SNP, EFFECT_ALLELE, NON_EFFECT_ALLELE, ZSCORE]
    if CHROMOSOME in d: keep_columns.append(CHROMOSOME)
    if POSITION in d: keep_columns.append(POSITION)
    if BETA in d: keep_columns.append(BETA)
    if SE in d: keep_columns.append(SE)
    if PVALUE in d: keep_columns.append(PVALUE)
    d = d[keep_columns]
    return d

def _rename_columns(d, gwas_format):
    # List of columns to try to rename
    cols = [ (COLUMN_SNP, SNP), (COLUMN_EFFECT_ALLELE, EFFECT_ALLELE), (COLUMN_NON_EFFECT_ALLELE, NON_EFFECT_ALLELE),
        (COLUMN_CHROMOSOME, CHROMOSOME), (COLUMN_POSITION, POSITION), (COLUMN_SE, SE),
        (COLUMN_BETA, BETA), (COLUMN_BETA_SIGN, BETA_SIGN), (COLUMN_OR, OR),
        (COLUMN_ZSCORE, ZSCORE), (COLUMN_PVALUE, PVALUE)]

    rename = {}
    for column, name in cols:
        if column in gwas_format:
            if gwas_format[column] in d:
                rename[gwas_format[column]] = name
            else:
                logging.info("Reading GWAS: Column %s not found", gwas_format[column])

    if len(rename):
        d = d.rename(columns=rename)

    return d

def _ensure_columns(d, input_pvalue_fix):
    if d.shape[0] == 0:
        if OR in d: d[BETA] = None
        if BETA_SIGN in d: d[BETA_SIGN] = None
        d[ZSCORE] = None
        return d

    d[EFFECT_ALLELE] = d[EFFECT_ALLELE].str.upper()
    d[NON_EFFECT_ALLELE] = d[NON_EFFECT_ALLELE].str.upper()

    if OR in d:
        logging.log(9, "Calculating beta from odds ratio")
        beta = _or_to_beta(d[OR])
        d[BETA] = beta

    if BETA_SIGN in d:
        b = d[BETA_SIGN]
        b = b.apply(lambda x: 1.0 if x == "+" else -1.0)
        d[BETA_SIGN] = b

    _ensure_z(d, input_pvalue_fix)

    d[ZSCORE] = numpy.array(d[ZSCORE], dtype=numpy.float32)
    return d

_numeric_columns = [BETA, OR, SE, PVALUE, ZSCORE]
def _enforce_numeric_columns(d):
    for column in _numeric_columns:
        if column in d:
            a = d[column]
            if a.dtype == numpy.object:
                a = [str(x) for x in a]
                a = [GWASSpecialHandling.sanitize_component(x) for x in a]
            d[column] = numpy.array(a, dtype=numpy.float64)
    return d

def _ensure_z(d, input_pvalue_fix):
    if ZSCORE in d:
        logging.log(9, "Using declared zscore")
        return d

    z = None

    if PVALUE in d:
        logging.log(9, "Calculating zscore from pvalue")
        z = _z_from_p(d, input_pvalue_fix)
    elif SE in d and BETA in d:
        logging.info("Calculating zscore from se and beta")
        z = d[BETA] / d[SE]

    if z is None: raise Exceptions.ReportableException("Couldn't get zscore from GWAS")
    d[ZSCORE] = z
    return d

def _z_from_p(d, input_pvalue_fix):
    p = d[PVALUE].values
    if numpy.any(p == 0):
        logging.warning("Encountered GWAS pvalues equal to zero. This might be caused by numerical resolution. Please consider using another scheme such as -beta- and -se- columns, or checking your input gwas for zeros.")

    s = _beta_sign(d)
    abs_z = -stats.norm.ppf(p / 2)

    if numpy.any(numpy.isinf(abs_z)) and input_pvalue_fix:
        logging.warning("Applying thresholding to divergent zscores. You can disable this behavior by using '--input_pvalue_fix 0' in the command line")
        the_min = numpy.min(p[numpy.logical_and(numpy.isfinite(abs_z),p != 0)])
        if input_pvalue_fix < the_min:
            the_min = input_pvalue_fix
        fix_z = -stats.norm.ppf(the_min / 2)
        logging.warning("Using %f to fill in divergent zscores", fix_z)
        abs_z[numpy.isinf(abs_z)] = fix_z

    z = abs_z * s
    return z

def _beta_sign(d):
    b = None
    if BETA in d:
        logging.log(9, "Acquiring sign from beta")
        b = numpy.sign(d[BETA])
    elif BETA_SIGN in d:
        logging.log(9, "Acquiring sign")
        b = d[BETA_SIGN]
        b = b.apply(lambda x: 1.0 if (x =="+" or x==1.0) else -1.0)
    if b is None: raise Exceptions.ReportableException("No beta sign in GWAS")
    return b

def _or_to_beta(odd):
    if numpy.any(numpy.where(odd < 0)):
        raise Exceptions.InvalidArguments("Odd Ratios include negative values.")
    if numpy.any(numpy.where(odd == 0)):
        logging.warning("Odd Ratios column holds some [0] values")
    return numpy.log(odd)

########################################################################################################################
# General purpose methods

def extract(gwas, snps):
    """Assumes that argument snps are there in the gwas."""
    s = gwas[SNP]
    index = [s[s == x].index[0] for x in snps]
    g = gwas.iloc[index]
    return g

def get_data_from_gwas(gwas, snp):
    if not snp in gwas[SNP]:
        return None
    g = gwas[gwas[SNP] == snp].iloc[0]
    return g

def process_gwas(path, callback, skip_header=True):
    with gzip.open(path) as file:
        for i,line in enumerate(file):
            if i==0 and skip_header: continue
            callback(line)

def get_header(gwas_path):
    with gzip.open(gwas_path) as file:
        header = file.readline().strip()
        return header

def get_snp_header_index(header, col_name):
    header = header.split()
    return header.index(col_name)

class GWASLineMappedCollector(object):
    def __init__(self, index_key=None):
        self.index_key = index_key
        self.collected = {}
        self.keys = []

    def __call__(self, line):
        comps = line.strip().split()
        k = comps[self.index_key]
        self.collected[k] = comps
        self.keys.append(k)
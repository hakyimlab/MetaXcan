import pandas #likely to go away for a streaming approach
import numpy
import logging
import os
import gzip
import scipy.stats as stats

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

#format of raw results
class GWASF(object):
    SNP=0
    CHROMOSOME=1
    POSITION=2
    NON_EFFECT_ALLELE=3
    EFFECT_ALLELE=4
    ZSCORE=5

def load_gwas(input_path, gwas_format, strict=True):
    """
    Attempts to read a GWAS summary statistics file, and load it into a uniform format,
    in a pandas dataframe.

    :param input_path: path to file containing GWAS summary statistics.
    :param gwas_format: dictionary specifying GWAS format column mapping.
        For example
    :return:
    """
    logging.info("Reading input gwas: %s", input_path)
    d = pandas.read_table(input_path)

    logging.info("Processing input gwas")
    d = _rename_columns(d, gwas_format)

    if not SNP in d:
        raise Exceptions.ReportableException("A valid SNP column name must be provided")

    #keep only rsids
    d = d[d[SNP].str.contains("rs")]

    if strict:
        d = _ensure_columns(d)
        d = _keep_gwas_columns(d)

    return d

def _or_to_beta(odd):
    return numpy.log(odd)

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

def _ensure_columns(d):
    if OR in d:
        logging.log(9, "Calculating beta from odds ratio")
        beta = _or_to_beta(d[OR])
        d[BETA] = beta

    if BETA_SIGN in d:
        b = d[BETA_SIGN]
        b = b.apply(lambda x: 1.0 if x == "+" else -1.0)
        d[BETA_SIGN] = b

    _ensure_z(d)

    d[ZSCORE] = numpy.array(d[ZSCORE], dtype=numpy.float32)
    return d

def _ensure_z(d):
    if ZSCORE in d:
        logging.log(9, "Using declared zscore")
        return

    z = None
    if PVALUE in d:
        logging.log(9, "Calculating zscore from pvalue")
        p = d[PVALUE]
        s = _beta_sign(d)
        abs_z = -stats.norm.ppf(p/2)
        z = abs_z*s
    elif SE in d and BETA in d:
        logging.info("Calculating zscore from se and beta")
        z= d[BETA]/d[SE]

    if z is None: raise Exceptions.ReportableException("Couldn't get zscore from GWAS")
    d[ZSCORE] = z
    return d

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
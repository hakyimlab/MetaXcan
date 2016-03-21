__author__ = 'heroico'

import gzip
import logging
import math
import scipy.stats as stats
import KeyedDataSet
import Utilities
import Exceptions
import os

class GWASTF(object):
    """GWAS file format"""
    SNP = 0
    A1 = 1
    A2 = 2
    FRQ = 3
    INFO = 4
    OR_BETA = 5
    SE = 6
    P = 7

    VALID_ALLELES =  ["A", "T", "C", "G"]

class GWASDosageFileIterator(object):
    """Exists mostly because gWAS output is not a CSV, and has variable whitespace between values"""
    def __init__(self, path=None, compressed=True, separator=None, callback=None):
        self.path = path
        self.callback = callback
        self.compressed = compressed
        self.separator = separator

    def iterateOverFile(self):
        class CallbackWrapper(object):
            def __init__(self, callback, separator):
                self.callback = callback
                self.separator = separator

            def __call__(self, i, line):
                comps = line.split(self.separator) if self.separator else line.split()
                line =  " ".join(comps)
                row = line.split()
                self.callback(row)

        file_iterator = Utilities.FileIterator(self.path, header="", compressed=self.compressed)
        callback = CallbackWrapper(self.callback, self.separator)
        file_iterator.iterate(callback)

class GWASFileFormat(object):
    def __init__(self, file_path, compressed, separator=None):
        if not os.path.isfile(file_path):
            raise Exceptions.BadFilename(file_path)
        header = None
        if compressed:
            with gzip.open(file_path, 'rb') as file:
                header = file.readline().strip()
        else:
            with open(file_path) as file:
                header = file.readline().strip()
        comps = header.split(separator) if separator else header.split()
        if "" in comps:
            comps.remove("")
        if len(comps) == 1: #didn't split the input string. Did we parse it with wrong arguments?
            raise Exception("Couldn't process input file. Please check for compressed status, or data field separator.")
        self.header_comps = [x for x in comps if x is not ""]
        self.SNP = None
        self.A1 = None
        self.A2 = None
        self.FRQ = None
        self.INFO = None
        self.OR = None
        self.BETA = None
        self.BETA_SIGN = None
        self.SE = None
        self.P = None
        self.BETA_ZSCORE = None

    def addSNPColumn(self, snp_column_name):
        if not snp_column_name in self.header_comps:
            raise NameError("SNP column name -%s- not found" % snp_column_name)

        self.SNP = self.header_comps.index(snp_column_name)

    def addSEColumn(self, se_column_name):
        if not se_column_name in self.header_comps:
            raise NameError("SE column name -%s- not found" % se_column_name)

        self.SE = self.header_comps.index(se_column_name)

    def addA1Column(self, A1_column_name):
        if not A1_column_name in self.header_comps:
            raise NameError("A1 column name -%s- not found" % A1_column_name)

        self.A1 = self.header_comps.index(A1_column_name)

    def addA2Column(self, A2_column_name):
        if not A2_column_name in self.header_comps:
            raise NameError("A1 column name -%s- not found" % A2_column_name)

        self.A2 = self.header_comps.index(A2_column_name)

    def addFrequencyColumn(self, frequency_column_name):
        if not frequency_column_name in self.header_comps:
            raise NameError("frequency column name -%s- not found" % frequency_column_name)

        self.FRQ = self.header_comps.index(frequency_column_name)

    def addORColumn(self, or_column_name):
        if not or_column_name in self.header_comps:
            raise NameError("OR column name -%s- not found" % or_column_name)

        self.OR = self.header_comps.index(or_column_name)

    def addBetaColumn(self, beta_column_name):
        if not beta_column_name in self.header_comps:
            raise NameError("beta column name -%s- not found" % beta_column_name)

        self.BETA = self.header_comps.index(beta_column_name)

    def addBetaSignColumn(self, beta_sign_column_name):
        if not beta_sign_column_name in self.header_comps:
            raise NameError("'beta sign' column name -%s- not found" % beta_sign_column_name)

        self.BETA_SIGN = self.header_comps.index(beta_sign_column_name)

    def addBetaZScoreColumn(self, beta_zscore_column_name):
        if not beta_zscore_column_name in self.header_comps:
            raise NameError("'beta zscore' column name -%s- not found" % beta_zscore_column_name)

        self.BETA_ZSCORE = self.header_comps.index(beta_zscore_column_name)

    def addPValueColumn(self, pvalue_column_name):
        if not pvalue_column_name in self.header_comps:
            raise NameError("'pvalue' column name '%s' not found", pvalue_column_name)

        self.P = self.header_comps.index(pvalue_column_name)

    @classmethod
    def fileFormatFromArgs(cls, file, args):
        file_format = GWASFileFormat(file, args.compressed, args.separator)
        if args.or_column:
            file_format.addORColumn(args.or_column)
        if args.beta_column:
            file_format.addBetaColumn(args.beta_column)
        if args.a1_column:
            file_format.addA1Column(args.a1_column)
        if args.a2_column:
            file_format.addA2Column(args.a2_column)
        if args.snp_column:
            file_format.addSNPColumn(args.snp_column)
        if args.frequency_column:
            file_format.addFrequencyColumn(args.frequency_column)
        if args.se_column:
            file_format.addSEColumn(args.se_column)
        if args.beta_zscore_column:
            file_format.addBetaZScoreColumn(args.beta_zscore_column)
        if args.beta_sign_column:
            file_format.addBetaSignColumn(args.beta_sign_column)
        if args.pvalue_column:
            file_format.addPValueColumn(args.pvalue_column)
        return file_format

class GWASSNPInfoLineCollector(object):
    A1=0
    A2=1
    OR_BETA=2
    FREQ=3

    def __init__(self, rsids=[], values=[]):
        self.rsids = rsids
        self.values = values
        self.beta = None
        self.ses = None
        self.beta_z = None
        self.sigma = None

    def __call__(self, row):
        a1 = row[GWASTF.A1].upper()
        a2 = row[GWASTF.A2].upper()
        OR_BETA = row[GWASTF.OR_BETA]
        freq=row[GWASTF.FRQ]
        self.rsids.append(row[GWASTF.SNP])
        self.values.append((a1,a2,OR_BETA, freq))

    def reset(self):
        self.rsids = []
        self.values = []

BETA = "beta"
BETA_SE = "beta_se"
BETA_SE_TO_Z = "beta_se_to_z"
Z = "z"
BETA_SIGN_P = "beta_sign_p"
BETA_P = "beta_p"

def _scheme(scheme, file_format):
    s = None
    if scheme == BETA:
        if bool(file_format.OR) == bool(file_format.BETA):
            raise Exception("Provide either 'beta' or 'or', but not both")
        s = _BETA_Scheme()
    elif scheme == BETA_SE:
        if bool(file_format.OR) == bool(file_format.BETA):
            raise Exception("Provide either 'beta' or 'or', but not both")
        if not file_format.SE:
            raise Exception("SE must be provided")
        s = _BETA_SE_Scheme()
    elif scheme == BETA_SE_TO_Z:
        if bool(file_format.OR) == bool(file_format.BETA):
            raise Exception("Provide either 'beta' or 'or', but not both")
        if not file_format.SE:
            raise Exception("SE must be provided")
        s = _BETA_SE_TO_Z_Scheme()
    elif scheme == Z:
        if not file_format.BETA_ZSCORE:
            raise Exception("'beta zscore column' must be provided")
        s = _BETA_Z_Scheme()
    elif scheme == BETA_P:
        if bool(file_format.OR) == bool(file_format.BETA):
            raise Exception("Provide either 'beta' or 'or', but not both")
        if not file_format.P:
            raise Exception("pvalue is required")
        s = _BETA_PVALUE_Scheme()
    elif scheme == BETA_SIGN_P:
        if not file_format.BETA_SIGN:
            raise Exception("'beta_sign_column' is required")
        if not file_format.P:
            raise Exception("pvalue is required")
        s = _BETA_SIGN_PVALUE_Scheme()
    else:
        raise Exception("Unknown scheme.")
    return s

class _GWASLineScheme(object):
    def __call__(self, collector, row, file_format):
        collector.rsids.append(row[file_format.SNP])

        sigma = "NA"
        if collector.file_format.FRQ:
            freq = row[file_format.FRQ]
            if freq != "NA":
                f = float(freq)
                sigma = str(math.sqrt(2 * f * (1-f)))
            collector.sigma.append(sigma)

def betaFromRow(file_format, row):
    beta = "NA"
    if file_format.BETA:
        beta = row[file_format.BETA]
    elif file_format.OR and not file_format.BETA:
        OR = row[file_format.OR]
        if OR != "NA":
            try:
                beta = math.log(float(OR))
                beta = str(beta)
            except Exception as e:
                logging.log(9, "Error at beta from row: %s", str(e))
                beta = "NA"
    return beta

def betaSignFromRow(file_format, row):
    sign = "NA"
    if file_format.BETA_SIGN:
        s = row[file_format.BETA_SIGN]
        if s == "+" or s == "-":
           sign = s
    return sign

class _BETA_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_Scheme, self).__call__(collector, row, file_format)

        beta = betaFromRow(file_format, row)
        collector.beta.append(beta)

class _BETA_SE_Scheme(_BETA_Scheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_SE_Scheme, self).__call__(collector, row, file_format)

        se = row[file_format.SE]
        collector.ses.append(se)

class _BETA_SE_TO_Z_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_SE_TO_Z_Scheme, self).__call__(collector, row, file_format)

        beta = betaFromRow(file_format, row)
        se = row[file_format.SE]
        beta_z = "NA"
        if beta != "NA" and "se" != "NA":
            beta_z = str(float(beta)/float(se))
        collector.beta_z.append(beta_z)

class _BETA_Z_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_Z_Scheme, self).__call__(collector, row, file_format)
        b_z = row[file_format.BETA_ZSCORE]
        collector.beta_z.append(b_z)

class _BETA_PVALUE_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_PVALUE_Scheme, self).__call__(collector, row, file_format)
        z = "NA"
        p = row[file_format.P]

        if p == ".":
            p = "NA"

        if p != "NA":
            p = float(p)
            abs_z = -stats.norm.ppf(p/2)
            beta = betaFromRow(file_format, row)
            if beta != "NA":
                try:
                    s = 1 if float(beta) >= 0 else -1
                    z = abs_z * s
                except Exception as e:
                    logging.log(9, "Error converting number: %s", str(e))
        collector.beta_z.append(z)

class _BETA_SIGN_PVALUE_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_SIGN_PVALUE_Scheme, self).__call__(collector, row, file_format)
        z = "NA"
        p = row[file_format.P]
        if p == ".":
            p = "NA"

        if p != "NA":
            p = float(p)
            abs_z = -stats.norm.ppf(p/2)
            sign = betaSignFromRow(file_format, row)
            if sign != "NA":
                s = 1.0 if sign == "+" else -1.0
                z = abs_z * s
        collector.beta_z.append(z)

class GWASBetaLineCollector(object):
    def __init__(self, file_format, scheme):
        self.file_format = file_format
        self.scheme = _scheme(scheme, file_format)
        self.reset()

    def __call__(self, row):
        self.scheme(self, row, self.file_format)

    def reset(self):
        self.rsids = []
        self.beta = [] if (self.file_format.BETA or self.file_format.OR) else None
        self.ses = [] if self.file_format.SE else None
        self.sigma = [] if ( self.file_format.FRQ ) else None #more coming soon
        self.OR = True if self.file_format.OR else None
        self.file_format = self.file_format
        self.beta_z = []

class GWASWeightDBFilteredBetaLineCollector(GWASBetaLineCollector):
    def __init__(self, file_format, scheme, weight_db_logic=None):
        super(GWASWeightDBFilteredBetaLineCollector, self).__init__(file_format, scheme)
        self.weight_db_logic = weight_db_logic

    def __call__(self, row):
        file_format = self.file_format
        rsid = row[file_format.SNP]
        if self.weight_db_logic:
            if not rsid in self.weight_db_logic.genes_for_an_rsid:
                logging.log(6, "%s not in weight db", rsid)
                return

            entry = self.weight_db_logic.anEntryWithRSID(rsid)
            a1 = row[file_format.A1].upper()
            a2 = row[file_format.A2].upper()
            if not a1 in GWASTF.VALID_ALLELES or \
                not a2 in GWASTF.VALID_ALLELES:
                logging.log(6,"invalid alleles %s %s", a1, a2)
                return

            if entry.ref_allele == a2 and entry.eff_allele == a1:
                logging.log(7, "alleles are flipped for rsid %s", rsid)
                if file_format.BETA:
                    try:
                        beta = betaFromRow(file_format, row)
                        b = -float(beta)
                        row[file_format.BETA] = str(b)
                    except Exception as e:
                        logging.log(9, "error flipping allele %s: %s", rsid, str(e))
                        row[file_format.BETA] = "NA"
                elif file_format.OR:
                    try:
                        beta = betaFromRow(file_format, row)
                        OR = math.exp(-float(beta))
                        row[file_format.OR] = str(OR)
                    except Exception as e:
                        logging.log(9, "error flipping allele %s: %s", rsid, str(e))
                        row[file_format.OR] = "NA"
                elif file_format.BETA_SIGN:
                    try:
                        sign = betaSignFromRow(file_format, row)
                        if sign == "+":
                            sign == "-"
                        elif sign == "-":
                            sign == "+"
                        elif sign == "NA":
                            pass
                        else:
                            raise Exception("wrong sign")
                        row[file_format.BETA_SIGN] = sign
                    except Exception as e:
                        logging.log(9, "error flipping allele %s: %s", rsid, str(e))
                        row[file_format.OR] = "NA"
            elif not entry.ref_allele == a1 or not entry.eff_allele == a2:
                logging.log(6, "%s alleles dont match:(%s, %s)(%s, %s)",rsid, entry.ref_allele, entry.eff_allele, a1, a2)
                return

        super(GWASWeightDBFilteredBetaLineCollector, self).__call__(row)


class GWASDosageFileLoader(object):
    def __init__(self, path, compressed=True, separator=None, callback=None, file_format=None, scheme=None):
        self.path = path
        self.callback = callback
        self.file_format = file_format
        self.compressed = compressed
        self.scheme = scheme
        self.separator = separator

    def load(self):
        callback = self.callback
        if not callback:
            logging.info("Default Beta callback")
            callback = GWASBetaLineCollector(self.path, self.file_format, self.scheme)

        file_iterator = GWASDosageFileIterator(self.path, self.compressed, self.separator, callback)
        file_iterator.iterateOverFile()

        results = []
        if callback.beta:
            beta = KeyedDataSet.KeyedDataSet(name="beta", data=callback.beta, keys=callback.rsids)
            results.append(beta)

        if callback.ses:
            se = KeyedDataSet.KeyedDataSet(name="se", data=callback.ses, keys=callback.rsids)
            results.append(se)

        if callback.beta_z:
            beta_z = KeyedDataSet.KeyedDataSet(name="beta_z", data=callback.beta_z, keys=callback.rsids)
            results.append(beta_z)

        if callback.sigma:
            freq = KeyedDataSet.KeyedDataSet(name="sigma_l", data=callback.sigma, keys=callback.rsids)
            results.append(freq)

        if type(callback) is GWASSNPInfoLineCollector:
            values = KeyedDataSet.KeyedDataSet(name="values", data=callback.values, keys=callback.rsids)
            results.append(values)

        callback.reset()
        return results

__author__ = 'heroico'

import gzip
import logging
import math
import scipy.stats as stats
import KeyedDataSet
import Utilities
import Exceptions
from Exceptions import ReportableException
import os

class GWASTF(object):
    """GWAS file format"""
    SNP = 0
    EFFECT_ALLELE = 1
    OTHER_ALLELE = 2
    FRQ = 3
    INFO = 4
    OR_BETA = 5
    SE = 6
    P = 7

    VALID_ALLELES =  ["A", "T", "C", "G"]

class GWASDosageFileIterator(object):
    """Exists mostly because gWAS output is not a CSV, and has variable whitespace between values"""
    def __init__(self, path=None, compressed=True, separator=None, callback=None, skip_until_header=None):
        self.path = path
        self.callback = callback
        self.compressed = compressed
        self.separator = separator
        self.skip_until_header = skip_until_header

    def iterateOverFile(self):
        class CallbackWrapper(object):
            def __init__(self, callback, separator):
                self.callback = callback
                self.separator = separator

            def __call__(self, i, line):
                comps = line.split(self.separator) if self.separator else line.split()
                line =  " ".join(comps)
                if len(line) == 1:
                    logging.log(8, "Found GWAS row with one component. Is your file ok?")
                    return
                row = line.split()
                self.callback(row)

        file_iterator = Utilities.FileIterator(self.path, header="", compressed=self.compressed) \
                            if not self.skip_until_header else \
                        Utilities.FileIterator(self.path, header=self.skip_until_header, compressed=self.compressed, ignore_until_header=True)

        callback = CallbackWrapper(self.callback, self.separator)
        file_iterator.iterate(callback)

class GWASFileFormat(object):
    def __init__(self, file_path, compressed, separator=None, skip_until_header=None):
        if not os.path.isfile(file_path):
            raise Exceptions.BadFilename(file_path)
        self.file_path = file_path
        header = None
        if compressed:
            with gzip.open(file_path, 'rb') as file:
                header = fileHeader(file, skip_until_header)
        else:
            with open(file_path) as file:
                header = fileHeader(file, skip_until_header)
        comps = header.split(separator) if separator else header.split()
        if "" in comps:
            comps.remove("")
        if len(comps) == 1: #didn't split the input string. Did we parse it with wrong arguments?
            raise Exception("Couldn't process input file. Please check for compressed status, or data field separator.")
        self.header_comps = [x for x in comps if x is not ""]
        self.SNP = None
        self.OTHER_ALLELE = None
        self.EFFECT_ALLELE = None
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
            raise ReportableException("SNP column name -%s- not found in file '%s'. Is the file compressed?" % (snp_column_name, self.file_path))

        self.SNP = self.header_comps.index(snp_column_name)

    def addSEColumn(self, se_column_name):
        if not se_column_name in self.header_comps:
            raise ReportableException("SE column name -%s- not found in file '%s'. Is the file compressed?" % (se_column_name, self.file_path))

        self.SE = self.header_comps.index(se_column_name)

    def addOtherAlleleColumn(self, oa_column_name):
        if not oa_column_name in self.header_comps:
            raise ReportableException("-other allele- column name -%s- not found in file '%s'. Is the file compressed?" % (oa_column_name, self.file_path))

        self.OTHER_ALLELE = self.header_comps.index(oa_column_name)

    def addEffectAlleleColumn(self, A2_column_name):
        if not A2_column_name in self.header_comps:
            raise ReportableException("-effect allele- column name -%s- not found in file '%s'. Is the file compressed?" % (A2_column_name, self.file_path))

        self.EFFECT_ALLELE = self.header_comps.index(A2_column_name)

    def addFrequencyColumn(self, frequency_column_name):
        if not frequency_column_name in self.header_comps:
            raise ReportableException("frequency column name -%s- not found in file '%s'. Is the file compressed?" % (frequency_column_name, self.file_path))

        self.FRQ = self.header_comps.index(frequency_column_name)

    def addORColumn(self, or_column_name):
        if not or_column_name in self.header_comps:
            raise ReportableException("OR column name -%s- not found in file '%s'. Is the file compressed?" % (or_column_name, self.file_path))

        self.OR = self.header_comps.index(or_column_name)

    def addBetaColumn(self, beta_column_name):
        if not beta_column_name in self.header_comps:
            raise ReportableException("beta column name -%s- not found in file '%s'. Is the file compressed?" % (beta_column_name, self.file_path))

        self.BETA = self.header_comps.index(beta_column_name)

    def addBetaSignColumn(self, beta_sign_column_name):
        if not beta_sign_column_name in self.header_comps:
            raise ReportableException("beta sign column name -%s- not found in file '%s'. Is the file compressed?" % (beta_sign_column_name, self.file_path))

        self.BETA_SIGN = self.header_comps.index(beta_sign_column_name)

    def addBetaZScoreColumn(self, beta_zscore_column_name):
        if not beta_zscore_column_name in self.header_comps:
            raise ReportableException("beta zscore column name -%s- not found in file '%s'. Is the file compressed?" % (beta_zscore_column_name, self.file_path))

        self.BETA_ZSCORE = self.header_comps.index(beta_zscore_column_name)

    def addPValueColumn(self, pvalue_column_name):
        if not pvalue_column_name in self.header_comps:
            raise ReportableException("pvalue column name -%s- not found in file '%s'. Is the file compressed?" % (pvalue_column_name, self.file_path))

        self.P = self.header_comps.index(pvalue_column_name)

    @classmethod
    def fileFormatFromArgs(cls, file, args):
        file_format = GWASFileFormat(file, args.compressed, args.separator, args.skip_until_header)
        if args.or_column:
            file_format.addORColumn(args.or_column)
        if args.beta_column:
            file_format.addBetaColumn(args.beta_column)
        if args.other_allele_column:
            file_format.addOtherAlleleColumn(args.other_allele_column)
        if args.effect_allele_column:
            file_format.addEffectAlleleColumn(args.effect_allele_column)
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

def fileHeader( file, skip_until_header):
    header = None
    if skip_until_header:
        for i,l in enumerate(file):
            l = l.strip()
            if skip_until_header in l:
                header = skip_until_header
                break
        if not header:
            raise  Exceptions.InvalidArguments("Wrong header lookup for GWAS files")
    else:
        header = file.readline().strip()
    return header

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
        collector.rsids.append(snpFromRow(file_format, row))
        if collector.gather_alleles:
            collector.ref_allele.append(row[file_format.A1])
            collector.eff_allele.append(row[file_format.A2])

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
        # Quick hack for files that miught be in non EN locale
        if "," in beta:
            beta = beta.replace(",",".")
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

def betaZscoreFromRow(file_format, row):
    b_z = "NA"
    if file_format.BETA_ZSCORE:
        b_z = row[file_format.BETA_ZSCORE]
    return b_z

def pFromFrow(file_format, row):
    p = row[file_format.P]

    if p == ".":
        p = "NA"

    # Quick hack for files that miught be in non EN locale
    if "," in p:
        p = p.replace(",",".")
    return p

def snpFromRow(file_format, row):
    snp = row[file_format.SNP]
    #Hacky heuristic for rs numbers
    if not "rs" in snp:
        #snps like "chr1:123456"
        if "chr" and ":" in snp:
            snp = "rs" + snp.split(":")[1]

    return snp

def betaZscoreFromRow(file_format, row):
    b_z = "NA"
    if file_format.BETA_ZSCORE:
        b_z = row[file_format.BETA_ZSCORE]
    return b_z

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
        collector.beta.append(beta)

class _BETA_Z_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_Z_Scheme, self).__call__(collector, row, file_format)
        b_z = betaZscoreFromRow(file_format, row)
        collector.beta_z.append(b_z)
        b = betaFromRow(file_format, row)
        collector.beta.append(b)

class _BETA_PVALUE_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_PVALUE_Scheme, self).__call__(collector, row, file_format)
        z = "NA"
        p = pFromFrow(file_format, row)
        beta = betaFromRow(file_format, row)

        if p != "NA":
            p = float(p)
            abs_z = -stats.norm.ppf(p/2)
            if beta != "NA":
                try:
                    s = 1 if float(beta) >= 0 else -1
                    z = abs_z * s
                except Exception as e:
                    logging.log(9, "Error converting number: %s", str(e))
        collector.beta_z.append(z)
        collector.beta.append(beta)

class _BETA_SIGN_PVALUE_Scheme(_GWASLineScheme):
    def __call__(self, collector, row, file_format):
        super(_BETA_SIGN_PVALUE_Scheme, self).__call__(collector, row, file_format)
        z = "NA"
        p = pFromFrow(file_format, row)

        if p != "NA":
            p = float(p)
            abs_z = -stats.norm.ppf(p/2)
            sign = betaSignFromRow(file_format, row)
            if sign != "NA":
                s = 1.0 if sign == "+" else -1.0
                z = abs_z * s
        collector.beta_z.append(z)
        #nothing we can do, if we got the sign of beta, then there is probably no beta present
        if collector.beta:
            b = betaFromRow(file_format, row)
            collector.beta.append(b)

class GWASBetaLineCollector(object):
    def __init__(self, file_format, scheme, gather_alleles=False):
        self.file_format = file_format
        self.scheme = _scheme(scheme, file_format)
        self.gather_alleles = gather_alleles
        self.reset()

    def __call__(self, row):
        self.scheme(self, row, self.file_format)

    def reset(self):
        self.rsids = []
        self.beta = [] if (self.file_format.BETA or self.file_format.OR or self.file_format.BETA_ZSCORE) else None
        self.ses = [] if self.file_format.SE else None
        self.sigma = [] if ( self.file_format.FRQ ) else None #more coming soon
        #self.OR = True if self.file_format.OR else None
        self.file_format = self.file_format
        self.beta_z = []
        self.ref_allele = [] if self.gather_alleles else None
        self.eff_allele = [] if self.gather_alleles else None

class GWASWeightDBFilteredBetaLineCollector(GWASBetaLineCollector):
    def __init__(self, file_format, scheme, weight_db_logic=None, gather_alleles=False):
        super(GWASWeightDBFilteredBetaLineCollector, self).__init__(file_format, scheme, gather_alleles)
        self.weight_db_logic = weight_db_logic

    def __call__(self, row):
        file_format = self.file_format
        rsid = snpFromRow(file_format, row)
        if self.weight_db_logic:
            if not rsid in self.weight_db_logic.genes_for_an_rsid:
                logging.log(6, "%s not in weight db", rsid)
                return

            other_allele = row[file_format.OTHER_ALLELE].upper()
            effect_allele = row[file_format.EFFECT_ALLELE].upper()
            if not other_allele in GWASTF.VALID_ALLELES or \
                not effect_allele in GWASTF.VALID_ALLELES:
                logging.log(6,"invalid alleles %s %s", effect_allele, other_allele)
                return

            # The following works but is inappropriate. All entries for a given SNP have the same ref allele.
            # but bear in mind that we are using any entry.
            entry = self.weight_db_logic.anEntryWithRSID(rsid)
            if entry.ref_allele == effect_allele and entry.eff_allele == other_allele:
                logging.log(7, "alleles are flipped for rsid %s", rsid)

                if file_format.BETA_ZSCORE:
                    try:
                        z = betaZscoreFromRow(file_format, row)
                        z = -float(z)
                        row[file_format.BETA_ZSCORE] = str(z)
                    except Exception as e:
                        logging.log(9, "error flipping allele %s: %s", rsid, str(e))
                        row[file_format.BETA_ZSCORE] = "NA"

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
                        row[file_format.BETA_SIGN] = "NA"
            elif not entry.ref_allele == other_allele or not entry.eff_allele == effect_allele:
                logging.log(6, "%s alleles dont match:(%s, %s)(%s, %s)",rsid, entry.ref_allele, entry.eff_allele, other_allele, effect_allele)
                return

        super(GWASWeightDBFilteredBetaLineCollector, self).__call__(row)

#convenience wrappers to output as we read
def loadGWASAndStream(input_path, output_path, compressed=True, separator=None, skip_until_header=None, callback=None, file_format=None, scheme=None):
    if not callback:
        logging.info("Default Beta callback")
        callback = GWASBetaLineCollector(file_format, scheme)

    class OutputWrapper(object):
        def __init__(self, collector, output_file):
            self.collector = collector
            self.output_file = output_file
            self.writeHeader()
            self.collected_anything = False

        def writeHeader(self):
            columns = ["rsid"]
            if self.collector.ref_allele is not None: columns.append("ref_allele")
            if self.collector.eff_allele is not None: columns.append("eff_allele")
            if self.collector.beta is not None: columns.append("beta")
            if self.collector.ses is not None: columns.append("beta_se")
            if self.collector.sigma is not None: columns.append("sigma")
            columns.append("beta_z")
            header = "%s\n" % (" ".join(columns))
            self.output_file.write(header)

        def __call__(self, row):
            self.collector(row)
            if len(self.collector.rsids):
                o = [self.collector.rsids[0]]
                if self.collector.ref_allele is not None: o.append(self.collector.ref_allele[0])
                if self.collector.eff_allele is not None: o.append(self.collector.eff_allele[0])
                if self.collector.beta is not None: o.append(str(self.collector.beta[0]))
                if self.collector.ses is not None: o.append(str(self.collector.ses[0]))
                if self.collector.sigma is not None: o.append(str(self.collector.sigma[0]))
                o.append(str(self.collector.beta_z[0]))
                line = "%s\n" % (" ".join(o))
                self.output_file.write(line)
                self.collected_anything = True
            self.collector.reset()

    def do_output(callback, output_file, input_path,  compressed=True, separator=None, skip_until_header=None):
        wrapper = OutputWrapper(callback, output_file)
        file_iterator = GWASDosageFileIterator(input_path, compressed, separator, wrapper, skip_until_header)
        file_iterator.iterateOverFile()
        if not wrapper.collected_anything:
            logging.info("No snps from the tissue model found in the GWAS file")

    if compressed:
        with gzip.open(output_path, "wb") as output_file:
            do_output(callback, output_file, input_path, compressed, separator, skip_until_header)
    else:
        with open(output_path, "w") as output_file:
            do_output(callback, output_file, input_path, compressed, separator, skip_until_header)

RSID="rsid"
BETA="beta"
BETA_SE="se"
BETA_Z="beta_z"
SIGMA_l="sigma_l"

class GWASDosageFileLoader(object):
    def __init__(self, path, compressed=True, separator=None, skip_until_header=None, callback=None, file_format=None, scheme=None):
        self.path = path
        self.callback = callback
        self.file_format = file_format
        self.compressed = compressed
        self.scheme = scheme
        self.separator = separator
        self.skip_until_header = skip_until_header

    def load(self):
        callback = self.callback
        if not callback:
            logging.info("Default Beta callback")
            callback = GWASBetaLineCollector(self.file_format, self.scheme)

        file_iterator = GWASDosageFileIterator(self.path, self.compressed, self.separator, callback, self.skip_until_header)
        file_iterator.iterateOverFile()

        if callback.rsids:
            results = {RSID:callback.rsids}
            column_order = [RSID]

        if callback.beta:
            results[BETA] = callback.beta
            column_order.append(BETA)

        if callback.ses:
            results[BETA_SE] = callback.ses
            column_order.append(BETA_SE)

        if callback.beta_z:
            results[BETA_Z] = callback.beta_z
            column_order.append(BETA_Z)

        if callback.sigma:
            results[SIGMA_l] = callback.sigma
            column_order.append(SIGMA_l)

        callback.reset()
        return results, column_order

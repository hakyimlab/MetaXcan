import gzip
import logging
import re
import io
import numpy
import pandas

from .. import Exceptions

def gwas_data_source(path, snps=None, snp_column_name=None, skip_until_header=None, separator=None, handle_empty_columns=False):
    s = {}
    def _ogz(p):
        return io.TextIOWrapper(gzip.open(p, "r"), newline="")
    _o = _ogz if ".gz" in path else open
    with _o(path) as file:
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
            duplicated = [k for k,v in header_count.items() if v>1]
            logging.info("The input GWAS has duplicated columns: %s, will only use the first one in each case", str(duplicated))

        if handle_empty_columns:
            split_r = re.compile(separator) if separator is not None else re.compile("\s")

        for i,line in enumerate(file):
            if handle_empty_columns:
                line = line.replace("\n", "")
                comps = split_r.split(line)
            else:
                comps = line.strip().split(separator)

            #Yeah, there are those kinds of files
            if not len(comps) == len(header_comps):
                logging.log(8, "Found line with less components than headers, line %i", i)
                continue

            if snps and not comps[index] in snps:
                continue

            # Load only the first column if in presence of duplicated columns. Yuck!
            sentinel=set()
            for i,c in enumerate(comps):
                comp = header_comps[i]
                if comp in sentinel: continue
                sentinel.add(comp)
                s[comp].append(c)

        for c in header_comps:
            s[c] = numpy.array(pandas.to_numeric(s[c], errors='ignore'))

    return s

non_en_number = re.compile("^[-\+]?[0-9]*,{1}[0-9]+([eE]{1}[-\+]?[0-9]+)?$")
def sanitize_component(c):
    if non_en_number.match(c): c = c.replace(",",".")
    if c == "": c = None
    elif c == "NA": c = None
    elif c == ".": c = None
    elif c == "\\N": c = None
    elif c == "": c= None

    return c

def to_numeric(d, column):
    if column in d:
        a = [sanitize_component(x) for x in d[column]]
        d[column] = numpy.array(a, dtype=numpy.float64)

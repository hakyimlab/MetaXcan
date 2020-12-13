__author__ = 'heroico'

import os
import io
import json
import re
import logging
import gzip
from . import Exceptions
import pandas
import numpy

VALID_ALLELES = ["A", "T", "C", "G"]

def hapName(name):
    return name + ".hap.gz"

def legendName(name):
    return name + ".legend.gz"

def dosageName(name):
    return name + ".dosage.gz"

def dosageNamesFromFolder(folder):
    names = contentsWithRegexpFromFolder(folder, ".*.dosage.gz")
    if not names:
        names = contentsWithRegexpFromFolder(folder, ".*.dos.gz")
    if not names:
        names = contentsWithRegexpFromFolder(folder, "chr.*.txt.gz")
    #if not names:
    #    raise Exceptions.ReportableException("Couldn't identify any dosages from {}".format(folder))
    return names

def hapNamesFromFolder(folder):
    names = namesWithPatternFromFolder(folder, ".hap.gz")
    return names

def legendNamesFromFolder(folder):
    names = namesWithPatternFromFolder(folder, ".legend.gz")
    return names

def namesWithPatternFromFolder(folder, pattern):
    contents = os.listdir(folder)
    names = []
    for content in contents:
        if pattern in content:
            name = content.split(pattern)[0]
            names.append(name)
    return names

def contentsWithPatternsFromFolder(folder, patterns):
    try:
        contents = os.listdir(folder)
        paths = []
        for content in contents:
            matches = True
            for pattern in patterns:
                if not pattern in content:
                    matches = False
                    break

            if matches:
                paths.append(content)
    except OSError:
        raise Exceptions.BadDirectory(folder)
    return paths

def contentsWithRegexpFromFolder(folder, regexp):
    if type(regexp) == str:
        regexp = re.compile(regexp)
    contents = os.listdir(folder)
    paths = [x for x in contents if regexp.match(x)] if regexp else contents
    return paths

def target_files(input_folder, file_filters=None):
    files = os.listdir(input_folder)
    if file_filters:
        patterns = [re.compile(x) for x in file_filters]
        files = [x for x in files if all([r.match(x) for r in patterns])]
    files = [os.path.join(input_folder, x) for x in files]
    return files

def samplesInputPath(path):
    samples_content = contentsWithPatternsFromFolder(path, [".sample"])
    if len(samples_content) == 0:
        samples_content = contentsWithPatternsFromFolder(path, ["samples.txt"])

    samples_path = None
    if len(samples_content) > 0:
        samples_file = samples_content[0]
        samples_path = os.path.join(path, samples_file)
    return  samples_path

def checkSubdirectorySanity(base, candidate):
    sane = True
    abs_base = os.path.realpath(os.path.abspath(base))
    abs_candidate = os.path.realpath(os.path.abspath(candidate))
    if abs_base == abs_candidate:
        return False

    rel_base_to_candidate = os.path.relpath(abs_base, abs_candidate)
    if not rel_base_to_candidate.startswith(os.pardir):
        return False

    return sane
open
class PercentReporter(object):
    def __init__(self, level, total, increment=10, pattern="%i percent complete"):
        self.level = level
        self.total = total
        self.increment = increment
        self.last_reported = 0
        self.pattern = pattern

    def update(self, i, text=None, force=False):
        percent = int(i*100.0/self.total)

        if force or  percent >= self.last_reported + self.increment:
            self.last_reported = percent

            if not text:
                text = self.pattern

            logging.log(self.level, text, percent)

def load_json(path):
    d = None
    with open(path) as file:
        d = json.load(file)
    return d

class FileIterator(object):
    def __init__(self, path, header=None, compressed = False, ignore_until_header = False):
        self.path = path
        self.compressed = compressed
        self.header = header
        self.ignore_until_header = ignore_until_header
        if ignore_until_header and not header:
            raise Exceptions.InvalidArguments("File iterator received conflicting header information")

    def iterate(self,callback=None):
        if self.compressed:
            with gzip.open(self.path, 'r') as file_object:
                self._iterateOverFile(io.TextIOWrapper(file_object, newline=""), callback)
        else:
            with open(self.path, 'r') as file_object:
                self._iterateOverFile(file_object, callback)

    def _iterateOverFile(self, file_object, callback):
        if self.ignore_until_header:
            self._ignore_until_header(file_object)
        else:
            if self.header is not None:
                line = file_object.readline().strip("\n")
                if len(self.header) and line != self.header:
                    raise Exceptions.MalformedInputFile(self.path, "Unexpected header")

        self._processFile(file_object, callback)

    def _ignore_until_header(self, file_object):
        if self.ignore_until_header and self.header:
            skip = True
            while skip:
                l = file_object.readline()
                if not l:
                    raise Exceptions.InvalidArguments("Wrong header lookup in %s" % (self.path,))
                l = l.strip()
                if self.header in l:
                    skip = False

    def _processFile(self, file_object, callback):
        if callback is not None:
            for i,line in enumerate(file_object):
                callback(i, line)

def open_any_plain_file(path):
    if "gz" in path:
        return io.TextIOWrapper(gzip.open(path), newline="")
    else:
        return open(path)

def generate_from_any_plain_file(path, skip_n=None):
    is_gz = ".gz" in path
    with open_any_plain_file(path) as f:
        if skip_n:
            for i in range(0, skip_n):
                f.readline()
        for l in f:
            yield l

import csv
class CSVFileIterator(FileIterator):
    def __init__(self, path, header=None, compressed = False, ignore_until_header = False, delimiter = " ", quotechar='"'):
        super(CSVFileIterator, self).__init__(path, header, compressed, ignore_until_header)
        self.delimiter = delimiter
        self.quotechar = quotechar

    def _processFile(self, file_object, callback):
        if callback is not None:
            reader = csv.reader(file_object, delimiter=self.delimiter, quotechar=self.quotechar)
            for i,row in enumerate(reader):
                callback(i, row)

def TS(string):
    """placeholder for string translation"""
    return string

def ensure_requisite_folders(path):
    folder = os.path.split(path)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)

def to_dataframe(data, columns,to_numeric=None, fill_na=None):
    data = list(zip(*data))
    if to_numeric:
        data = [pandas.to_numeric(x, errors=to_numeric) for x in data]
    if len(data) == 0: data = [[] for i in range(0, len(columns))]
    data = {columns[i]:data[i] for i in range(0, len(columns))}
    data = pandas.DataFrame(data)
    data = data[columns]
    if fill_na:
        data = data.fillna(fill_na)
    return data

def save_dataframe(d, path, mode="w", header=True):
    compression = "gzip" if "gz" in path else None
    ensure_requisite_folders(path)
    d.to_csv(path, header=header, mode=mode, compression=compression, sep="\t", index=False, na_rep="NA")

def save_table(data, path, mode="w", header=None):
    compression = "gzip" if "gz" in path else None
    def _ogz(p):
        return  io.TextIOWrapper(gzip.open(p, mode), newline="")
    _o = _ogz if ".gz" in path else open
    ensure_requisite_folders(path)

    def _to_line(comps):
        line = "\t".join([str(x) for x in comps]) + "\n"
        return line

    with _o(path, mode) as file:
        if header:
            file.write(_to_line(header))

        for d in data:
            line = _to_line(d)
            file.write(line)

def sub_batch(d, sub_batches, sub_batch):
    batches = numpy.array_split(range(0, d.shape[0]), sub_batches)
    batch = batches[sub_batch]
    d = d.iloc[batch]
    d = d.reset_index(drop=True)
    return d
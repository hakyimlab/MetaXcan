__author__ = 'heroico'

import os
import json
import Exceptions

VALID_ALLELES = ["A", "T", "C", "G"]

def hapName(name):
    return name + ".hap.gz"

def legendName(name):
    return name + ".legend.gz"

def dosageName(name):
    return name + ".dosage.gz"

def dosageNamesFromFolder(folder):
    names = namesWithPatternFromFolder(folder, ".dosage.gz")
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
    contents = os.listdir(folder)
    paths = [x for x in contents if regexp.match(x)] if regexp else contents
    return paths

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

import logging
class PercentReporter(object):
    def __init__(self, level, total, increment=10, pattern="%i percent complete"):
        self.level = level
        self.total = total
        self.increment = increment
        self.last_reported = 0
        self.pattern = pattern

    def update(self, i, text=None):
        percent = int(i*100.0/self.total)
        if percent == self.last_reported + self.increment:
            self.last_reported = percent

            if not text:
                text = self.pattern

            logging.log(self.level, text, percent)

def load_json(path):
    d = None
    with open(path) as file:
        d = json.load(file)
    return d

import gzip
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
            with gzip.open(self.path, 'rb') as file_object:
                self._iterateOverFile(file_object, callback)
        else:
            with open(self.path, 'rb') as file_object:
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
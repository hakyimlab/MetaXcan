__author__ = 'heroico'

import gzip
import io
import logging
from .DataSet import DataSet

class KeyedDataSet(DataSet):
    def __init__(self, name=None, index=None, data=[], keys=[]):
        super(KeyedDataSet, self).__init__(name, index, data)
        self.keys = keys
        self.values_by_key = {}
        if len(keys):
            self.values_by_key = {keys[i]:data[i] for i in range(0, len(keys))}

EMPTY = KeyedDataSet()

class KeyedDataSetFileUtilities(object):
#SAVE
    @classmethod
    def saveToFile(cls, file_path, keyed_data_set, key_header, value_header):
        with open(file_path, "w") as file:
            cls.writeContents(file, keyed_data_set, key_header, value_header)

    @classmethod
    def saveToCompressedFile(cls, file_path, keyed_data_set, key_header, value_header):
        with io.TextIOWrapper(gzip.open(file_path, 'w'), newline="") as file:
            cls.writeContents(file, keyed_data_set, key_header, value_header)

    @classmethod
    def writeContents(self, file, keyed_data_set, key_header, value_header):
        if key_header and value_header:
            file.write(key_header + " " + value_header + "\n")
        for i in range(0, len(keyed_data_set.keys)):
            key = keyed_data_set.keys[i]
            value = keyed_data_set.data[i]
            file.write(key + " " + str(value) + "\n")

#SAVE SEVERAL DATA SETS AT ONCE
    @classmethod
    def saveSetsToFile(cls, file_path, sets, key_header):
        with open(file_path, "w") as file:
            cls.writeSetsContent(file, sets, key_header)

    @classmethod
    def saveSetsToCompressedFile(cls, file_path, sets, key_header):
        with io.TextIOWrapper(gzip.open(file_path, 'w'), newline="") as file:
            cls.writeSetsContent(file, sets, key_header)

    @classmethod
    def writeSetsContent(self, file, sets, key_header):
        sets_count = len(sets)
        if key_header:
            names = [x.name for x in sets]
            values_header = " ".join(names)
            file.write(key_header + " " + values_header + "\n")

        if not len(sets):
            logging.info("No KeyedDataSets to save, skipping")
            return

        keys = sets[0].keys
        reference = set(keys)
        for i in range(1, sets_count):
            assert reference == set(sets[i].keys)

        for key in keys:
            values = [str(x.values_by_key[key]) for x in sets]
            line = " ".join(values)
            file.write(key + " " + line + "\n")

#LOAD
    @classmethod
    def loadFromFile(cls, file_path, sep= None, col = 1, header=None):
        with open(file_path, 'r') as file:
            keyed_data_set = cls.loadContents(file, file_path, sep, col, header)
            return keyed_data_set

    @classmethod
    def loadFromCompressedFile(cls, file_path, sep= None, col = 1, header=None):
        with io.TextIOWrapper(gzip.open(file_path, 'r'), newline="") as file:
            keyed_data_set = cls.loadContents(file, file_path, sep, col, header)
            return keyed_data_set

    @classmethod
    def loadContents(cls, file, file_path, separator=None, col=1, header=None):
        if header is not None:
            file_header = file.readline().strip()
            if len(header):
                if file_header != header:
                    logging.info("header issue:\n%s\n%s", file_header, header)
                assert file_header == header

        keys = []
        values = []
        for line in file:
            comps = line.strip().split(separator) if separator else line.strip().split()
            keys.append(comps[0])
            values.append(comps[col])

        keyed_data_set = KeyedDataSet(file_path, None, values, keys)
        return keyed_data_set

#LOAD SEVERAL DATA SETS FROM A FILE AT ONCE
    @classmethod
    def loadDataSetsFromFile(cls, file_path, sep= None, cols = [], header=None):
        with open(file_path, 'r') as file:
            sets = cls.loadDataSetsContent(file, sep, cols, header)
            return sets

    @classmethod
    def loadDataSetsFromCompressedFile(cls, file_path, sep= None, cols = [], header=None):
        with io.TextIOWrapper(gzip.open(file_path, 'r'), newline="") as file:
            sets = cls.loadDataSetsContent(file, sep, cols, header)
            return sets

    @classmethod
    def loadDataSetsContent(cls, file, separator=None, cols=[], header=None):
        names = None
        if header is not None:
            line = file.readline().strip("\n")
            if len(header):
                assert header == line
            names = line.split(separator) if separator else line.split()

        if len(cols) == 0:
            for i in range(0, len(names)):
                cols.append(i)

        value_sets = []
        for i,col in enumerate(cols):
            value_sets.append([])

        keys = []
        for line in file:
            comps = line.strip().split(separator) if separator else line.strip().split()
            keys.append(comps[0])
            for i,col in enumerate(cols):
                value_set = value_sets[i]
                value_set.append(comps[col])

        keyed_data_sets = []
        for i, col in enumerate(cols):
            value_set = value_sets[i]
            name = names[col] if names else None
            keyed_data_set = KeyedDataSet(name=name, index=col, data=value_set, keys=keys)
            keyed_data_sets.append(keyed_data_set)
        return keyed_data_sets

def setWithName(sets, name):
    result = None
    for candidate in sets:
        if candidate.name == name:
            result = candidate
            break
    return result

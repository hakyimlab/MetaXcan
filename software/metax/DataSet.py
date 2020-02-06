__author__ = 'heroico'

import gzip
import io
#
class DataSet(object):
    """A list of values for a thing we care about"""
    def __init__(self,name=None,index=None, data = []):
        self.data = data
        self.name = name
        self.index = index


class DataSetFileUtilities(object):
    @classmethod
    def loadFromFile(cls, data_file_name = None, header_name=None):
        data_set = None
        with open(data_file_name, 'r') as file:
            data_set = cls._loadDataSetFromFile(file, header_name)
            data_set.name = data_file_name
        return data_set

    @classmethod
    def loadFromCompressedFile(cls, data_file_name = None, header_name=None):
        data_set = None
        with io.TextIOWrapper(gzip.open(data_file_name, 'r'), newline="") as file:
            data_set = cls._loadDataSetFromFile(file, header_name)
            data_set.name = data_file_name
        return data_set

    @classmethod
    def _loadDataSetFromFile(cls,file, header_name):
        read_first_line = False
        data = []
        if header_name is not None:
            line = file.readline().strip("\n")
            assert line.strip("\n") == header_name

        for i,line in enumerate(file):
            data.append(line.strip("\n"))
        data_set = DataSet()
        data_set.data = data
        return data_set

#
class DataSetCollection(object):
    """A list of values for a thing we care about"""
    def __init__(self):
        self.sets = []


__author__ = 'heroico'


import sqlite3
import numpy
from .. import KeyedDataSet
from .. import Exceptions
import os

class DBLoaders(object):
    @classmethod
    def loadKeyedDataSetFromDB(cls, db_path, table_name,  key_column, value_column):
        if not os.path.exists(db_path):
            raise Exceptions.BadFilename(db_path)
        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()
        query = "SELECT "+key_column+", "+value_column+" FROM " + table_name
        try:
            entries = cursor.execute(query)
        except sqlite3.Error as e:
            raise Exceptions.InvalidDbFormat(db_path, e.args[0])

        rsids = []
        variances = []
        for entry in entries:
            rsids.append(entry[0])
            variances.append(entry[1])

        keyed_data_set = KeyedDataSet.KeyedDataSet(db_path, None, variances, rsids)
        return keyed_data_set

    @classmethod
    def loadVariancesFromDB(cls, db_path):
        keyed_data_set = cls.loadKeyedDataSetFromDB(db_path, "variances", "rsid", "var")
        return keyed_data_set

    @classmethod
    def loadCovarianceMatrix(cls, db_path, keys):
        if not os.path.exists(db_path):
            raise Exceptions.BadFilename(db_path)
        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()

        key_filter = {}
        valid_keys = []         # list of valid rsid 1s

        values = {}             # Lookup for all observed pairwise covariances
                                # rsx=>{rsy=>X,},
        def get_row(dict, key):
            row = None
            if key in dict:
                row = dict[key]
            else:
                row = {}
                dict[key] = row
            return row

        cursor.execute("SELECT rsid1, rsid2, covariance FROM covariances")
        results = cursor.fetchall()
        for result in results:
            rsid1 = result[0]
            rsid2 = result[1]
            value = result[2]

            if value == "NA":
                continue

            if not rsid1 in key_filter:
                key_filter[rsid1] = True
                valid_keys.append(rsid1)

            value = float(value)

            row_1 = get_row(values, rsid1)
            row_1[rsid2] = value

            row_2 = get_row(values, rsid2)
            row_2[rsid1] = value
        connection.close()

        valid_rows = []
        for i in range(0, len(valid_keys)):
            valid_row = []
            valid_rows.append(valid_row)
            for j in range(0, len(valid_keys)):
                key_i = valid_keys[i]
                key_j = valid_keys[j]

                value = values[key_i][key_j]
                valid_row.append(value)

        covariance_matrix = numpy.array(valid_rows)
        return covariance_matrix, valid_keys
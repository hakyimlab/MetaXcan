__author__ = 'heroico'


import sqlite3
import numpy
import KeyedDataSet

class DBLoaders(object):
    @classmethod
    def loadKeyedDataSetFromDB(cls, db_path, table_name,  key_column, value_column):
        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()
        query = "SELECT "+key_column+", "+value_column+" FROM " + table_name
        entries = cursor.execute(query)

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
        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()

        key_filter = {}
        valid_keys = []

        values = {}
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
        for i in xrange(0, len(valid_keys)):
            valid_row = []
            valid_rows.append(valid_row)
            for j in xrange(0, len(valid_keys)):
                key_i = valid_keys[i]
                key_j = valid_keys[j]

                value = values[key_i][key_j]
                valid_row.append(value)

        covariance_matrix = numpy.array(valid_rows)
        return covariance_matrix, valid_keys
__author__ = 'heroico'

import numpy
import logging
import Utilities
import Exceptions

class MTF(object):
    """Matrix table format"""
    GENE = 0
    RSID_1 = 1
    RSID_2 = 2
    VALUE = 3

def loadMatrixFromFile(file):
    class MatrixBuilder(object):
        def __init__(self):
            self.build_up = {}
            self.entries = {}
            self.latest_gene = None
            self.pending = []

        def __call__(self, i, row):
            gene = row[MTF.GENE]
            rsid1 = row[MTF.RSID_1]
            rsid2 = row[MTF.RSID_2]
            value = row[MTF.VALUE]

            if self.latest_gene != gene:
                self.processPending()
                self.latest_gene = gene

            if value == "NA":
                return

            self.pending.append((gene, rsid1, rsid2, value))

        def processPending(self):
            if not len(self.pending):
                return

            the_gene = self.pending[0][0]
            if the_gene in self.entries:
                raise Exceptions.MalformedInputFile(file,
                        "Snp Covariance Entries for genes must be contiguous but %s was found in two different places in the file." % (the_gene))

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

            for gene, rsid1, rsid2, value in self.pending:
                if not rsid1 in key_filter:
                    key_filter[rsid1] = True
                    valid_keys.append(rsid1)

                value = float(value)

                row_1 = get_row(values, rsid1)
                row_1[rsid2] = value

                row_2 = get_row(values, rsid2)
                row_2[rsid1] = value

            valid_rows = []
            for i in xrange(0, len(valid_keys)):
                valid_row = []
                valid_rows.append(valid_row)
                for j in xrange(0, len(valid_keys)):
                    key_i = valid_keys[i]
                    key_j = valid_keys[j]

                    value = values[key_i][key_j]
                    valid_row.append(value)

            self.pending = [] #pending data is not needed anymore
            covariance_matrix = numpy.array(valid_rows)
            self.entries[the_gene] = (covariance_matrix, valid_keys)

    builder = MatrixBuilder()
    loader = Utilities.CSVFileIterator(file, header="GENE RSID1 RSID2 VALUE", compressed=True)
    loader.iterate(builder)
    builder.processPending() #Builder read last entries, but since the gene didn't change, it didn't realize it has one last matrix to build
    return builder.entries

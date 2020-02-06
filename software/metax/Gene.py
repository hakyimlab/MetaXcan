__author__ = 'heroico'

import logging
from . import Utilities

class Gene(object):
    def __init__(self, name, base_position, chromosome_name, ens_id):
        self.name = name
        self.base_position = base_position
        self.chromosome_name = chromosome_name
        self.ens_id = ens_id

    @classmethod
    def loadFromDigest(cls, path, header=None, compressed=False, autosome_only=True):
        class DFT(object):
            COL_CHR=0
            COL_SIGN=1
            COL_BASE_POSITION=2
            COL_END_POSITION=3
            COL_ENSEMBLE_ID=4
            COL_NAME=5
            COL_WHATEVER=6
            COL_STATUS=7

        class GeneCollectorCallback(object):
            def __init__(self, autosome_only):
                self.autosome_only = autosome_only
                self.collected = {}
                self.collected_by_name = {}

            def __call__(self, i, row):
                chr = row[DFT.COL_CHR]
                chr_blah = chr.split("chr")[1]
                if self.autosome_only and chr_blah.lower() in ["x", "y"]:
                    return

                base_position = row[DFT.COL_BASE_POSITION]
                name = row[DFT.COL_NAME]

                ens_id = row[DFT.COL_ENSEMBLE_ID]
                if ens_id in self.collected:
                    logging.info("Gene %s already there", name)
                    return

                gene = Gene(name, base_position, chr_blah, ens_id)
                self.collected[ens_id] = gene

                if name in self.collected_by_name:
                    logging.info("Gene name %s already there", name)
                    return
                self.collected_by_name[name] = gene

        file_iterator = Utilities.CSVFileIterator(path, header, compressed, delimiter="\t")
        callback = GeneCollectorCallback(autosome_only)
        file_iterator.iterate(callback)
        return callback.collected, callback.collected_by_name

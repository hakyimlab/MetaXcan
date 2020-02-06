__author__ = 'heroico'

from . import DataSet

class DataSetSNP(DataSet.DataSet):
    """" A dat aset form an snp """
    def __init__(self, name=None, index=None, data=[], position=None, ref_allele=None, eff_allele=None):
        super(DataSetSNP, self).__init__(name, index, data)
        self.position = position
        self.ref_allele = ref_allele
        self.eff_allele = eff_allele

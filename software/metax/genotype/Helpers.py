__author__ = "alvaro barbeira"

from .. import Exceptions
from ..misc import KeyedDataSource

def gtex_snp(gtex_snp_file, gtex_release_version):
    if not gtex_release_version:
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "VariantID", "RS_ID_dbSNP142_CHG37p13", numeric=False)
    elif gtex_release_version.lower() == "v7":
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "variant_id", "rs_id_dbSNP147_GRCh37p13", numeric=False)
    elif gtex_release_version.lower() == "v8":
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "variant_id", "rs_id_dbSNP150_GRCh38p7", numeric=False)
    elif gtex_release_version.lower() == "model_training_v7":
        gtex_snp = KeyedDataSource.load_data(gtex_snp_file, "varID", "rsid_dbSNP150", numeric=False)
    else:
        raise Exceptions.InvalidArguments("Unsupported GTEx release version")
    return gtex_snp
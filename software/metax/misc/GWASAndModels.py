import logging
import os
import re
import pandas

from .. import Constants
from .. import PredictionModel
from ..gwas import Utilities as GWASUtilities

def align_data_to_alleles(data, base, left_on, right_on):
    EA, NEA = Constants.EFFECT_ALLELE, Constants.NON_EFFECT_ALLELE
    EA_BASE, NEA_BASE = EA+"_BASE", NEA+"_BASE"
    merged = pandas.merge(data, base, left_on=left_on, right_on=right_on, suffixes=("", "_BASE"))

    alleles_1 = pandas.Series([set(e) for e in zip(merged[EA], merged[NEA])])
    alleles_2 = pandas.Series([set(e) for e in zip(merged[EA_BASE], merged[NEA_BASE])])
    eq = alleles_1 == alleles_2
    merged = merged[eq]
    if eq.shape[0] == 0:
        return merged

    flipped = merged[EA] != merged[EA_BASE]
    Z = Constants.ZSCORE
    if Z in merged:
        merged.loc[flipped, Z] = - merged.loc[flipped, Z]
    B = Constants.BETA
    if B in merged:
        merged.loc[flipped, B] = - merged.loc[flipped, B]

    merged.loc[flipped, EA] = merged.loc[flipped, EA_BASE]
    merged.loc[flipped, NEA] = merged.loc[flipped, NEA_BASE]

    return merged

def gwas_model_intersection(args):
    gwas= GWASUtilities.load_plain_gwas_from_args(args)
    paths = PredictionModel._model_paths(args.models_folder, args.models_name_filter)
    PF = PredictionModel.WDBQF
    intersection = set()
    for db_path in sorted(paths):
        logging.log(9, "loading %s", db_path)
        model = PredictionModel.load_model(db_path)
        base = model.weights[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE]].drop_duplicates()
        b = align_data_to_alleles(gwas, base, Constants.SNP, PF.K_RSID)
        intersection.update(b[Constants.SNP])
    return intersection



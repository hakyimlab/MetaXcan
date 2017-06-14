import logging
import os
import re
import pandas

from .. import Utilities
from .. import Constants
from .. import PredictionModel
from ..gwas import GWAS
from ..gwas import Utilities as GWASUtilities

def align_data_to_alleles(data, base, left_on, right_on):
    EA, NEA = Constants.EFFECT_ALLELE, Constants.NON_EFFECT_ALLELE
    EA_BASE, NEA_BASE = EA+"_BASE", NEA+"_BASE"
    merged = pandas.merge(data, base, left_on=left_on, right_on=right_on, suffixes=("", "_BASE"))

    alleles_1 = pandas.Series([set(e) for e in zip(merged[EA], merged[NEA])])
    alleles_2 = pandas.Series([set(e) for e in zip(merged[EA_BASE], merged[NEA_BASE])])
    eq = alleles_1 == alleles_2
    merged = merged[eq]

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

def load_gwas(args):
    regexp = re.compile(args.gwas_file_pattern) if args.gwas_file_pattern else  None
    gwas_format = GWASUtilities.gwas_format_from_args(args)
    GWAS.validate_format_basic(gwas_format)
    GWAS.validate_format_for_strict(gwas_format)

    names = Utilities.contentsWithRegexpFromFolder(args.gwas_folder, regexp)
    names.sort()  # cosmetic, because different filesystems/OS yield folders in different order
    load_from = [os.path.join(args.gwas_folder,x) for x in names]
    sep = '\s+' if args.separator is None else args.separator
    if args.skip_until_header:
        _l = lambda x: GWASUtilities.gwas_filtered_source(x, skip_until_header=args.skip_until_header, separator=args.separator)
        load_from = [_l(x) for x in load_from]
    files = [GWAS.load_gwas(x, gwas_format, sep=sep) for x in load_from]
    gwas = pandas.concat(files)
    return gwas

def gwas_model_intersection(args):
    gwas= load_gwas(args)
    paths = PredictionModel._model_paths(args.models_folder)
    PF = PredictionModel.WDBQF
    intersection = set()
    for db_path in sorted(paths):
        logging.log(9, "loading %s", db_path)
        model = PredictionModel.load_model(db_path)
        base = model.weights[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE]].drop_duplicates()
        b = align_data_to_alleles(gwas, base, Constants.SNP, PF.K_RSID)
        intersection.update(b[Constants.SNP])
    return intersection



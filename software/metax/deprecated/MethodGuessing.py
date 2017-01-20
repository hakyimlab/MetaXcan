__author__ = 'heroico'

import os
import logging
import gzip
import ZScoreCalculation
import GWASUtilities
import Normalization
import Exceptions


def chooseZscoreSchemeFromFiles(folder, beta_contents, covariance_entries, weight_db_logic):
    if len(beta_contents) == 0:
        raise Exceptions.ReportableException("No snp's beta data found. Please check your beta files and/or command line arguments.")

    beta_content = beta_contents[0]
    beta_path = os.path.join(folder, beta_content)
    zscore_scheme = None
    normalization_scheme = None
    with gzip.open(beta_path) as content:
        header = content.readline().strip()
        # So. If beta_z is present, just go for "modified formula" with reference variance.
        # Any other option has to be specifically chosen.
        if "beta_z" in header:
            zscore_scheme = ZScoreCalculation.BETA_Z_SIGMA_REF
            normalization_scheme = Normalization.NONE
        elif "beta" in header and "sigma_l" in header:
            zscore_scheme = ZScoreCalculation.METAXCAN
            normalization_scheme = _chooseNormalization(header)
        elif "beta" in header and not "sigma_l" in header:
            zscore_scheme = ZScoreCalculation.METAXCAN_FROM_REFERENCE
            normalization_scheme = _chooseNormalization(header)
        else:
            raise Exception("Couldn't infer data from beta file header")
    logging.info("Chose zscore scheme '%s' and normalization '%s'", zscore_scheme, normalization_scheme)
    zscore_calculation = ZScoreCalculation. ZScoreScheme(zscore_scheme)
    normalization = Normalization.normalizationScheme(normalization_scheme, covariance_entries, weight_db_logic)
    return zscore_calculation, normalization

def _chooseNormalization(header):
    normalization_scheme = Normalization.NONE
    if "se" in header and "sigma_l" in header:
        normalization_scheme = Normalization.FROM_PHENO
    elif "se" in header and not "sigma_l" in header:
        normalization_scheme = Normalization.FROM_REFERENCE
    return normalization_scheme

def chooseGWASProcessingScheme(args, input_path):
    if args.scheme is not None:
        return args.scheme

    # just guess among the PValue schemes.
    # The other schemes are available only through specific selection, if you know what you are doing.
    scheme = None
    if args.pvalue_column and (args.beta_column or args.or_column):
        scheme = GWASUtilities.BETA_P
    elif args.pvalue_column and args.beta_sign_column:
        scheme = GWASUtilities.BETA_SIGN_P
    elif args.beta_zscore_column:
        scheme = GWASUtilities.Z
    else:
        raise Exception("Couldn't infer data processing choice. PValue and one from [beta, OR, sign of beta] is needed. ZScore of snp against phenotype is a viable alternative.")

    logging.info("Selected scheme '%s' for %s", scheme, input_path)
    return scheme

def chooseGWASCallback(file_format, scheme, weight_db_logic):
    if weight_db_logic:
        callback = GWASUtilities.GWASWeightDBFilteredBetaLineCollector(file_format, scheme, weight_db_logic=weight_db_logic)
    else:
        callback = GWASUtilities.GWASBetaLineCollector(file_format, scheme, gather_alleles=True)
    return callback
"""Choosing which of a model db's variant id columns to match the input on.

A PredictDB model db carries more than one usable variant id column: MASHR and
GTEx v8 dbs have both `rsid` and `varID` (`chr9_137032730_G_A_b38`), older ones
only `rsid`. `PredictionModel.ModelDB` matches on exactly one of them, chosen
by --model_db_snp_key, and picking the wrong one is silent: the merge succeeds
with nothing in it and the run finishes with an empty or all-NA results file.

The same is true one layer over, for the covariance: PredictDB ships MASHR
covariances keyed by varID, so a MASHR model loaded by rsid matches the GWAS
and then finds nothing in the covariance.

This module samples all three sides, works out which column can actually be
matched, and says out loud what it picked. It only ever picks a column that is
already in the db -- no id is constructed or rewritten here.
"""

import logging
import os
import re
import sqlite3

import pandas

from . import SnpOverlapDiagnostics
from .. import Utilities

# Columns of a weights table that are never variant ids.
_NOT_AN_ID = ("gene", "weight", "ref_allele", "eff_allele")

DEFAULT_SNP_KEY = "rsid"

UNRECOGNIZED = SnpOverlapDiagnostics.UNRECOGNIZED


########################################################################################################################
# Sampling the three sides. Each returns None when it can't tell, and the
# caller then leaves everything alone: guessing from a failed read is worse
# than the mismatch we are trying to catch.

def model_db_id_columns(path):
    """The weights-table columns of a model db that could serve as a matching key."""
    if not path or not os.path.exists(path):
        return []
    try:
        connection = sqlite3.connect(path)
        try:
            columns = [x[1] for x in connection.execute("PRAGMA table_info(weights);")]
        finally:
            connection.close()
    except sqlite3.Error as e:
        logging.log(9, "Could not inspect model db columns: %s", str(e))
        return []
    return [x for x in columns if x.lower() not in _NOT_AN_ID]

def sample_model_ids(path, column, n=1000):
    """A sample of one weights-table column's values, to tell its id format."""
    try:
        connection = sqlite3.connect(path)
        try:
            # the column name comes from the db's own schema, never from user input
            q = 'SELECT "{}" FROM weights LIMIT {};'.format(column, int(n))
            return [x[0] for x in connection.execute(q) if x[0] is not None]
        finally:
            connection.close()
    except sqlite3.Error as e:
        logging.log(9, "Could not sample model db column %s: %s", column, str(e))
        return []

def sample_gwas_ids(path, snp_column, separator=None, n=1000):
    """A sample of the GWAS's variant ids, read straight off the file.

    Deliberately reads only the id column and only the first n rows: this runs
    before the GWAS is loaded for real, and must not cost anything noticeable.
    """
    if separator is None or separator == "ANY_WHITESPACE":
        separator = '\s+'
    try:
        d = pandas.read_table(path, sep=separator, usecols=[snp_column], nrows=n)
    except Exception as e:
        logging.log(9, "Could not sample GWAS variant ids from %s: %s", path, str(e))
        return None
    # x == x drops NaN
    return [x for x in d[snp_column].values if x is not None and x == x]

def sample_covariance_ids(path, id_column="RSID1", n=1000):
    """A sample of the covariance's variant ids, read straight off the file."""
    if not path or not os.path.exists(path):
        return None
    try:
        d = pandas.read_table(path, sep='\s+', usecols=[id_column], nrows=n)
    except Exception as e:
        logging.log(9, "Could not sample covariance variant ids from %s: %s", path, str(e))
        return None
    return [x for x in d[id_column].values if x is not None and x == x]


########################################################################################################################
# The decision itself, kept free of io so it can be tested on formats alone.

_WIKI = "See https://github.com/hakyimlab/MetaXcan/wiki for guidance."

def choose_model_snp_key(model_formats, current_key, gwas_format, covariance_format=None):
    """Which of the model db's id columns to match the GWAS on.

    model_formats maps each candidate column of the weights table to the id
    format its values are in. Returns (key, message): key is None when the
    current key should stand, message is None when there is nothing worth
    telling the user.

    Only a column whose format matches the GWAS's is ever chosen, and only when
    the user did not name a key explicitly (that is the caller's business).
    An unrecognized format on any side means we can't reason about it, so we
    say nothing rather than guess.
    """
    if gwas_format == UNRECOGNIZED or current_key not in model_formats:
        return None, None

    current_format = model_formats[current_key]
    cov_known = covariance_format not in (None, UNRECOGNIZED)

    if current_format == gwas_format:
        if cov_known and covariance_format != current_format:
            return None, _unfixable_conflict_message(current_key, current_format, gwas_format, covariance_format)
        return None, None

    candidates = sorted(k for k, f in model_formats.items() if f == gwas_format)
    if not candidates:
        return None, None

    if cov_known and covariance_format != gwas_format:
        return None, _unfixable_conflict_message(current_key, current_format, gwas_format, covariance_format)

    return candidates[0], _switch_message(current_key, current_format, candidates[0], gwas_format)

def _switch_message(current_key, current_format, new_key, gwas_format):
    return (
        "The model db's -{}- column holds {} ids but the GWAS's variant ids look like {}, "
        "so matching on -{}- would have found almost nothing. Matching on the db's -{}- column "
        "instead, which is in the same format as the GWAS.\n"
        "Pass --model_db_snp_key {} to make this explicit, or --model_db_snp_key {} to force the "
        "original behavior."
    ).format(current_key, current_format, gwas_format, current_key, new_key, new_key, current_key)

def _unfixable_conflict_message(current_key, current_format, gwas_format, covariance_format):
    return (
        "The GWAS's variant ids look like {} and the covariance's look like {}, so no single "
        "column of the model db can match both: matching the GWAS leaves the covariance empty "
        "(every gene computed from zero SNPs, all-NA results) and matching the covariance leaves "
        "the GWAS empty. Currently matching on -{}-, which holds {} ids.\n"
        "Things to check:\n"
        "- Use a GWAS column in the same format as the covariance if your input has one "
        "(harmonized GWAS usually carry both an rsID and a chr_pos_ref_alt id).\n"
        "- Make sure the covariance you passed to --covariance is the one distributed with this "
        "model.\n"
    ).format(gwas_format, covariance_format, current_key, current_format) + _WIKI

def should_keep_non_rsid(gwas_format):
    """Whether the GWAS loader's rsid-only filter has to be lifted.

    GWAS.load_gwas drops every variant whose id doesn't contain "rs" unless
    --keep_non_rsid is passed, so a chr_pos_ref_alt-keyed GWAS loses all of its
    rows before anything else gets a chance to go wrong.
    """
    return gwas_format not in (UNRECOGNIZED, SnpOverlapDiagnostics.RSID)


########################################################################################################################
# Glue: same decision, driven off a parsed command line.

def _gwas_path(args):
    if getattr(args, "gwas_file", None):
        return args.gwas_file
    folder = getattr(args, "gwas_folder", None)
    if not folder:
        return None
    pattern = getattr(args, "gwas_file_pattern", None)
    names = Utilities.contentsWithRegexpFromFolder(folder, re.compile(pattern) if pattern else None)
    if not names:
        return None
    names.sort()
    return os.path.join(folder, names[0])

def resolve_snp_key_arguments(args):
    """Line up --model_db_snp_key and --keep_non_rsid with the input's actual ids.

    Mutates args in place and logs what it changed. Returns the list of
    messages, for tests. Anything the user asked for explicitly is left alone:
    an explicit --model_db_snp_key is never overridden, only diagnosed.
    """
    messages = []

    model_db_path = getattr(args, "model_db_path", None)
    if not model_db_path:
        return messages

    if getattr(args, "skip_until_header", None):
        # the header isn't where a plain read expects it, and reproducing
        # GWASSpecialHandling's parsing just to peek isn't worth it
        logging.log(9, "Skipping snp key resolution: --skip_until_header in use")
        return messages

    path = _gwas_path(args)
    snp_column = getattr(args, "snp_column", None)
    if not path or not snp_column:
        return messages

    gwas_ids = sample_gwas_ids(path, snp_column, getattr(args, "separator", None))
    if not gwas_ids:
        return messages
    gwas_format = SnpOverlapDiagnostics.dominant_id_format(gwas_ids)

    covariance_ids = sample_covariance_ids(getattr(args, "covariance", None))
    covariance_format = SnpOverlapDiagnostics.dominant_id_format(covariance_ids) if covariance_ids else None

    user_key = getattr(args, "model_db_snp_key", None)
    current_key = user_key if user_key else DEFAULT_SNP_KEY
    columns = model_db_id_columns(model_db_path)
    model_formats = {c: SnpOverlapDiagnostics.dominant_id_format(sample_model_ids(model_db_path, c)) for c in columns}

    key, message = choose_model_snp_key(model_formats, current_key, gwas_format, covariance_format)
    if key is None:
        if message:
            messages.append(message)
    elif not user_key:
        args.model_db_snp_key = key
        messages.append(message)
    else:
        # the user named a key; report the disagreement but do as told
        messages.append(
            "Using --model_db_snp_key {} as requested, but the GWAS's variant ids look like {}, "
            "which matches the model db's -{}- column instead. Drop --model_db_snp_key to let it "
            "be chosen from the data.".format(user_key, gwas_format, key))

    if should_keep_non_rsid(gwas_format) and not getattr(args, "keep_non_rsid", False):
        args.keep_non_rsid = True
        messages.append(
            "The GWAS's variant ids look like {}, not rsIDs, so every one of them would have been "
            "dropped on load. Enabling --keep_non_rsid.".format(gwas_format))

    for m in messages:
        logging.warning(m)
    return messages

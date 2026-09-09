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
from .. import PredictionModel
from .. import Utilities

# Columns of a weights table that are never variant ids.
_NOT_AN_ID = ("gene", "weight", "ref_allele", "eff_allele")

DEFAULT_SNP_KEY = "rsid"

UNRECOGNIZED = SnpOverlapDiagnostics.UNRECOGNIZED


########################################################################################################################
# Sampling the three sides. Each returns None when it can't tell, and the
# caller then leaves everything alone: guessing from a failed read is worse
# than the mismatch we are trying to catch.

def _query_model_db(path, query):
    """Runs a read-only query against a model db, or None if it can't be read.

    Goes through os.path.exists first because sqlite3.connect happily creates
    an empty db at a path that doesn't exist, which would turn a mistyped
    --model_db_path into a silent empty model.
    """
    if not path or not os.path.exists(path):
        return None
    try:
        connection = sqlite3.connect(path)
        try:
            return list(connection.execute(query))
        finally:
            connection.close()
    except sqlite3.Error as e:
        logging.log(9, "Could not read model db %s: %s", path, str(e))
        return None

def model_db_id_columns(path):
    """The weights-table columns of a model db that could serve as a matching key."""
    rows = _query_model_db(path, "PRAGMA table_info(weights);")
    if rows is None:
        return []
    return [x[1] for x in rows if x[1].lower() not in _NOT_AN_ID]

def sample_model_ids(path, column, n=1000):
    """A sample of one weights-table column's values, to tell its id format."""
    # the column name comes from the db's own schema, never from user input
    rows = _query_model_db(path, 'SELECT "{}" FROM weights LIMIT {};'.format(column, int(n)))
    return [x[0] for x in rows if x[0] is not None] if rows else []

def model_id_translation(path, from_column, to_column):
    """A variant id lookup built from the model db's own weights table.

    Both columns describe the same row, hence the same variant with the same
    alleles, so this renames ids and nothing else -- no position is parsed and
    no allele is inferred. Ids that name more than one variant in the target
    column (the same rsID at a multi-allelic site, say) are dropped rather than
    resolved arbitrarily.
    """
    if from_column == to_column:
        return {}
    # column names come from the db's own schema, never from user input
    rows = _query_model_db(path, 'SELECT DISTINCT "{}", "{}" FROM weights;'.format(from_column, to_column))
    if not rows:
        return {}

    pairs = [x for x in rows if x[0] is not None and x[1] is not None]
    targets = PredictionModel.patch_variant_names([x[1] for x in pairs])

    translation, ambiguous = {}, set()
    for (source, _), target in zip(pairs, targets):
        if source in translation and translation[source] != target:
            ambiguous.add(source)
        translation[source] = target
    for source in ambiguous:
        del translation[source]
    if ambiguous:
        logging.log(9, "Dropped %d ambiguous variant ids from the translation", len(ambiguous))
    return translation

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

def choose_keys(model_formats, current_key, gwas_format, covariance_format=None, allow_key_switch=True):
    """Which db column to key the model on, and which one the GWAS's ids are in.

    model_formats maps each candidate column of the weights table to the id
    format its values hold. Returns (model_key, gwas_key, message):

    - model_key is the column to match on, or None to leave the current one be.
      It is chosen to agree with the covariance when the covariance's format is
      known, because a model the covariance can't be looked up in produces an
      all-NA results file no matter how well the GWAS matched.
    - gwas_key is the column the GWAS's own ids are in, when that isn't the
      model key -- the GWAS then has to be translated from that column to the
      key, which the db itself can do since it carries both.
    - message is None when there is nothing worth telling the user.

    An unrecognized format on any side means we can't reason about it, so we
    say nothing rather than guess.
    """
    if gwas_format == UNRECOGNIZED or current_key not in model_formats:
        return None, None, None

    # first column holding each format, so the choice is deterministic
    by_format = {}
    for column, fmt in sorted(model_formats.items()):
        by_format.setdefault(fmt, column)

    gwas_column = by_format.get(gwas_format)
    if gwas_column is None:
        # nothing in the db is in the GWAS's format; neither a column choice nor
        # a translation through the db can help, and the overlap diagnostics
        # will say so with the actual numbers
        return None, None, None

    cov_known = covariance_format not in (None, UNRECOGNIZED)
    if cov_known:
        target = by_format.get(covariance_format)
        if target is None:
            return None, None, _foreign_covariance_message(model_formats, covariance_format)
    else:
        target = gwas_column

    if not allow_key_switch and target != current_key:
        message = _explicit_key_message(current_key, target, gwas_format, covariance_format if cov_known else None)
        target = current_key
    else:
        message = None

    model_key = target if target != current_key else None
    gwas_key = gwas_column if gwas_column != target else None

    if message is None and (model_key or gwas_key):
        message = _plan_message(current_key, model_formats[current_key], target, gwas_column,
                                gwas_format, covariance_format if cov_known else None)
    return model_key, gwas_key, message

def _plan_message(current_key, current_format, target, gwas_column, gwas_format, covariance_format):
    if gwas_column == target:
        return (
            "The model db's -{}- column holds {} ids but the GWAS's variant ids look like {}, "
            "so matching on -{}- would have found almost nothing. Matching on the db's -{}- column "
            "instead, which is in the same format as the GWAS.\n"
            "Pass --model_db_snp_key {} to make this explicit, or --model_db_snp_key {} to force "
            "the original behavior."
        ).format(current_key, current_format, gwas_format, current_key, target, target, current_key)

    return (
        "The GWAS's variant ids look like {} but the covariance is keyed by {} ids, so the model "
        "has to be matched on its -{}- column for the covariance to be found at all -- otherwise "
        "every gene is computed from zero SNPs and every result is NA. Translating the GWAS's ids "
        "from the db's -{}- column to -{}-, which is a lookup in the model db's own weights table; "
        "alleles are matched afterwards as usual.\n"
        "Pass --model_db_snp_key {} --gwas_snp_key {} to make this explicit."
    ).format(gwas_format, covariance_format, target, gwas_column, target, target, gwas_column)

def _explicit_key_message(current_key, target, gwas_format, covariance_format):
    return (
        "Using --model_db_snp_key {} as requested, but the GWAS's ids look like {}{} and the db's "
        "-{}- column is the one that lines up. Drop --model_db_snp_key to let it be chosen from "
        "the data."
    ).format(current_key, gwas_format,
             ", the covariance's like {}".format(covariance_format) if covariance_format else "",
             target)

def _foreign_covariance_message(model_formats, covariance_format):
    return (
        "The covariance's variant ids look like {}, which is not the format of any id column in "
        "the model db ({}). Every gene will be computed from zero SNPs and the results will be "
        "all NA.\n"
        "- Check that the covariance you passed to --covariance is the one distributed with this "
        "model.\n"
    ).format(covariance_format,
             ", ".join("{}: {}".format(k, v) for k, v in sorted(model_formats.items()))) + _WIKI

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
    """Line up --model_db_snp_key, --gwas_snp_key and --keep_non_rsid with the input's actual ids.

    Mutates args in place and logs what it changed. Returns the list of
    messages, for tests. Anything the user asked for explicitly is left alone:
    an explicit --model_db_snp_key is never overridden, only diagnosed, and an
    explicit --gwas_snp_key or --snp_map_file means the ids are already being
    mapped by hand and none of this applies.
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

    if getattr(args, "snp_map_file", None) or getattr(args, "gwas_snp_key", None):
        logging.log(9, "Skipping snp key resolution: variant ids are being mapped explicitly")
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

    model_key, gwas_key, message = choose_keys(model_formats, current_key, gwas_format,
                                               covariance_format, allow_key_switch=not user_key)
    if message:
        messages.append(message)
    if model_key:
        args.model_db_snp_key = model_key
    if gwas_key:
        args.gwas_snp_key = gwas_key

    if should_keep_non_rsid(gwas_format) and not getattr(args, "keep_non_rsid", False):
        args.keep_non_rsid = True
        messages.append(
            "The GWAS's variant ids look like {}, not rsIDs, so every one of them would have been "
            "dropped on load. Enabling --keep_non_rsid.".format(gwas_format))

    for m in messages:
        logging.warning(m)
    return messages

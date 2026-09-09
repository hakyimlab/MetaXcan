import itertools
import re

RSID = "rsid"
VARID = "chr_pos_ref_alt"
CHR_POS = "chr:pos"
UNRECOGNIZED = "unrecognized"

_FORMATS = (
    (RSID, re.compile(r"^rs\d+$", re.IGNORECASE)),
    (VARID, re.compile(r"^(chr)?[0-9XYM]+(T)?_\d+_[ACGTN]+_[ACGTN]+(_b\d+)?$", re.IGNORECASE)),
    (CHR_POS, re.compile(r"^(chr)?[0-9XYM]+(T)?[:_]\d+([:_][ACGTN]+[:_][ACGTN]+)?$", re.IGNORECASE)),
)

def classify_snp_id(snp_id):
    for name, pattern in _FORMATS:
        if pattern.match(str(snp_id)):
            return name
    return UNRECOGNIZED

def dominant_id_format(snps, sample_n=1000):
    """The most common variant id format among an arbitrary sample of snps."""
    counts = {}
    for snp in itertools.islice(iter(snps), sample_n):
        k = classify_snp_id(snp)
        counts[k] = counts.get(k, 0) + 1
    if not counts:
        return UNRECOGNIZED
    return max(sorted(counts), key=lambda k: counts[k])

def example_ids(snps, n=5, prefer_format=None):
    """A few ids to show the user, favoring ones in the format we claim is dominant.

    Without the preference a plain sort leads with whichever minority format
    sorts first, which undercuts a message about the ids looking like something
    else.
    """
    ids = sorted(snps)
    if prefer_format is None:
        return ids[:n]
    preferred = list(itertools.islice((x for x in ids if classify_snp_id(x) == prefer_format), n))
    return preferred or ids[:n]

def formats_disagree(format_1, format_2):
    """Whether two id formats are recognized and different.

    Two unrecognized formats are not called a mismatch: they may well be the
    same in-house scheme, and a false alarm here is worse than a missed one.
    """
    if UNRECOGNIZED in (format_1, format_2):
        return False
    return format_1 != format_2

def check_snp_overlap_from_counts(model_snps, matched_count, min_pct=1.0, sample_n=5):
    """Compares a model's SNP set against a count of matched variants.

    Used where the "other side" (e.g. streamed genotype input) is never
    materialized as a full set, so only a running count of matches is
    available. Returns (overlap_pct, message); message is None unless
    overlap_pct is below min_pct.
    """
    total = len(model_snps)
    if total == 0:
        return 0.0, None

    pct = 100.0 * matched_count / total
    if pct >= min_pct:
        return pct, None

    model_sample = sorted(model_snps)[:sample_n]
    message = (
        "Only {:.2f}% of the model's SNPs were found in the genotype input "
        "({} of {}). This usually means the model and genotype data use "
        "different variant id formats or genome builds, and the predicted "
        "expression will be mostly empty.\n"
        "Example model SNP ids: {}\n"
        "Things to check:\n"
        "- If your genotype variant ids don't match the model's rsids "
        "directly, use --variant_mapping (or --on_the_fly_mapping) to map "
        "genotype ids to model rsids.\n"
        "- If ids match but alleles don't (--skip_palindromic drops "
        "ambiguous variants), check your allele coding.\n"
        "- Check that your genotype data and model use the same genome "
        "build (GTEx v8 models are hg38); use --liftover otherwise.\n"
        "- If your genotype data deliberately covers only part of the genome "
        "(a single chromosome, say), low overlap is expected and you can "
        "ignore this.\n"
        "- As a rule of thumb, overlap below ~80% is suspect; single digits "
        "to teens usually means a mismatch. See "
        "https://github.com/hakyimlab/MetaXcan/wiki for guidance."
    ).format(pct, matched_count, total, model_sample)
    return pct, message

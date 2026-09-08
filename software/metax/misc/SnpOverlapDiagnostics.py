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
        "- As a rule of thumb, overlap below ~80% is suspect; single digits "
        "to teens usually means a mismatch. See "
        "https://github.com/hakyimlab/MetaXcan/wiki for guidance."
    ).format(pct, matched_count, total, model_sample)
    return pct, message

import logging
import pandas
import os
import numpy
from scipy import stats

from .. import Constants
from .. import Utilities
from .. import MatrixManager
from ..PredictionModel import WDBQF, WDBEQF, load_model, dataframe_from_weight_data
from ..misc import DataFrameStreamer
from ..misc import SnpOverlapDiagnostics
from . import AssociationCalculation

class SimpleContext(AssociationCalculation.Context):
    def __init__(self, gwas, model, covariance):
        self.gwas = gwas
        self.model = model
        self.covariance = covariance

    def get_weights(self, gene):
        w = self.model.weights
        w = w[w.gene == gene]
        return w

    def get_covariance(self, gene, snps):
        return self.covariance.get(gene, snps, strict_whitelist=False)


    def get_n_in_covariance(self, gene):
        return self.covariance.n_ids(gene)

    def get_gwas(self, snps):
        g = self.gwas
        g = g[g[Constants.SNP].isin(snps)]
        return g

    def get_model_snps(self):
        return set(self.model.weights.rsid)

    def get_gwas_snps(self):
        return set(self.gwas[Constants.SNP])

    def get_covariance_snps(self):
        return self.covariance.snps()

    def get_data_intersection(self):
        return _data_intersection(self.model, self.gwas)

    def provide_calculation(self, gene):
        w = self.get_weights(gene)
        gwas = self.get_gwas(w[WDBQF.K_RSID].values)

        i = pandas.merge(w, gwas, left_on="rsid", right_on="snp")
        if not Constants.BETA in i: i[Constants.BETA] = None
        i = i[[Constants.SNP, WDBQF.K_WEIGHT, Constants.ZSCORE, Constants.BETA]]

        snps, cov = self.get_covariance(gene, i[Constants.SNP].values)

        # fast subsetting and aligning
        d_columns = i.columns.values
        if snps is not None and len(snps):
            d = {x[0]: x for x in i.values}
            d = [d[snp] for snp in snps]
            d = list(zip(*d))
            d = {d_columns[i]:d[i] for i in range(0, len(d_columns))}
            i = pandas.DataFrame(d)
        else:
            i = pandas.DataFrame(columns=d_columns)
        return len(w.weight), i, cov, snps

    def get_model_info(self):
        return self.model.extra

class OptimizedContext(SimpleContext):
    def __init__(self, gwas, model, covariance, MAX_R):
        self.covariance = covariance
        self.genes, self.weight_data, self.snps_in_model = _prepare_weight_data(model, MAX_R)
        self.gwas_data = _prepare_gwas_data(gwas)
        self.extra = model.extra
        self.last_gene = None
        self.data_cache = None
        self.pedantic = MAX_R is None

    def _get_weights(self, gene):
        w = self.weight_data[gene]
        w = {x[WDBQF.RSID]:x[WDBQF.WEIGHT] for x in w}
        return w

    def get_weights(self, gene):
        w = self.weight_data[gene]
        w = dataframe_from_weight_data(list(zip(*w)))
        return w

    def get_model_snps(self):
        return set(self.snps_in_model)

    def get_gwas_snps(self):
        return set(self.gwas_data.keys())

    def _get_gwas(self, snps):
        snps = set(snps)
        g = self.gwas_data
        g = [g[x] for x in snps if x in g]
        g = {x[0]:(x[1], x[2]) for x in g}
        return g

    def get_gwas(self, snps):
        snps = set(snps)
        g = self.gwas_data
        g = [g[x] for x in snps if x in g]
        if len(g):
            g = list(zip(*g))
            g = pandas.DataFrame({Constants.SNP:g[0], Constants.ZSCORE:g[1], Constants.BETA:g[2]})
        else:
            g = pandas.DataFrame(columns=[Constants.SNP, Constants.ZSCORE, Constants.BETA])
        return g

    def get_data_intersection(self):
        return _data_intersection_3(self.weight_data, self.gwas_data, self.extra.gene.values, self.pedantic)

    def provide_calculation(self, gene):
        if gene != self.last_gene:
            #dummy while(True) to emulate go/to
            while True:
                w = self._get_weights(gene)
                gwas = self._get_gwas(list(w.keys()))
                type = [str, numpy.float64, numpy.float64, numpy.float64]
                columns = [Constants.SNP, WDBQF.K_WEIGHT, Constants.ZSCORE, Constants.BETA]
                d = {x: v for x, v in w.items() if x in gwas}

                snps, cov = self.get_covariance(gene, list(d.keys()))
                if snps is None:
                    d = pandas.DataFrame(columns=columns)
                    self.data_cache = len(w), d, cov, snps
                    self.last_gene = gene
                    break

                d = [(x, w[x], gwas[x][0], gwas[x][1]) for x in snps]
                d = list(zip(*d))
                if len(d):
                    d = {columns[i]:numpy.array(d[i], dtype=type[i]) for i in range(0,len(columns))}
                else:
                    d = {columns[i]:numpy.array([]) for i in range(0,len(columns))}

                self.data_cache = len(w), d, cov, snps
                self.last_gene = gene
                break

        return self.data_cache

    def get_model_info(self):
        return self.extra


def _data_intersection(model, gwas):
    weights = model.weights
    k = pandas.merge(weights, gwas, how='inner', left_on="rsid", right_on="snp")
    genes = k.gene.drop_duplicates().values
    snps = k.rsid.drop_duplicates().values
    return genes, snps

def _data_intersection_2(weight_data, gwas_data):
    genes = set()
    snps = set()
    for gene, entries in weight_data.items():
        gs = list(zip(*entries))[WDBQF.RSID]
        for s in gs:
            if s in gwas_data:
                genes.add(gene)
                snps.add(s)
    return genes, snps

def _data_intersection_3(weight_data, gwas_data, gene_list, pedantic):
    genes = list()
    _genes = set()
    snps =set()
    for gene in gene_list:
        if not gene in weight_data:
            if pedantic:
                logging.warning("Issues processing gene %s, skipped", gene)
            continue
        gs = list(zip(*weight_data[gene]))[WDBQF.RSID]
        for s in gs:
            if s in gwas_data:
                if not gene in _genes:
                    _genes.add(gene)
                    genes.append(gene)
                snps.add(s)

    return genes, snps

_ID_ADVICE = (
    "- If your model uses MASHR (varID-keyed) variants, pass "
    "--model_db_snp_key varID and make sure --snp_column points at a "
    "matching chr_pos_ref_alt(_build) column, or provide --snp_map_file.\n"
    "- If your GWAS uses non-rsID variant ids, pass --keep_non_rsid, which "
    "otherwise drops every id that doesn't look like an rsID.\n"
)

_BUILD_ADVICE = (
    "- Check that your GWAS and model use the same genome build "
    "(GTEx v8 models are hg38; liftover otherwise) and compatible "
    "reference panels.\n"
)

_WIKI_ADVICE = "See https://github.com/hakyimlab/MetaXcan/wiki for guidance."

def check_snp_overlap(model_snps, gwas_snps, min_pct=1.0, sample_n=5):
    """Compares the model's SNP set against the GWAS's SNP set.

    Returns (overlap_pct, message). message is None unless the two look
    mismatched, in which case it is a human-readable diagnostic meant to help
    identify SNP id/build mismatches (e.g. rsid vs varID) that would otherwise
    silently produce a near-empty result set.

    A low overlap percentage on its own is not enough to call a mismatch: a
    GWAS covering a single chromosome legitimately overlaps a genome-wide
    model by a few percent. So the primary signal is the variant id format
    differing between the two sides, which stays true regardless of how much
    of the genome the GWAS covers. min_pct is kept as a floor for the case
    formats can't see, notably an hg19/hg38 mismatch where both sides are
    chr_pos_ref_alt but no position ever lines up.
    """
    total = len(model_snps)
    if total == 0:
        return 0.0, None

    overlap = model_snps & gwas_snps
    pct = 100.0 * len(overlap) / total

    model_format = SnpOverlapDiagnostics.dominant_id_format(model_snps)
    gwas_format = SnpOverlapDiagnostics.dominant_id_format(gwas_snps)
    format_mismatch = SnpOverlapDiagnostics.formats_disagree(model_format, gwas_format)

    if pct >= min_pct and not format_mismatch:
        return pct, None

    model_examples = SnpOverlapDiagnostics.example_ids(model_snps, sample_n, model_format)
    examples = (
        "Example model SNP ids: {}\n"
        "Example GWAS SNP ids: {}\n"
    ).format(model_examples, SnpOverlapDiagnostics.example_ids(gwas_snps, sample_n, gwas_format))

    if not gwas_snps:
        message = (
            "No GWAS variants at all survived loading, so none of the model's {} "
            "SNPs could be matched and the results file will be empty.\n"
            "Example model SNP ids: {}\n"
            "Things to check:\n"
        ).format(total, model_examples) + _ID_ADVICE + _BUILD_ADVICE + _WIKI_ADVICE
    elif format_mismatch:
        message = (
            "The model's variant ids and the GWAS's variant ids are in different "
            "formats (model looks like {}, GWAS looks like {}), so they can't be "
            "matched: only {:.2f}% of the model's SNPs were found in the GWAS data "
            "({} of {}) and the results file will be mostly empty.\n"
        ).format(model_format, gwas_format, pct, len(overlap), total) + examples + \
            "Things to check:\n" + _ID_ADVICE + _WIKI_ADVICE
    else:
        message = (
            "Only {:.2f}% of the model's SNPs were found in the GWAS data ({} of {}), "
            "so most genes will have no result. The variant ids are in the same format "
            "on both sides ({}), so this is not an id format mismatch.\n"
        ).format(pct, len(overlap), total, model_format) + examples + (
            "Things to check:\n"
            "- If your GWAS deliberately covers only part of the genome (a single "
            "chromosome, say), low overlap is expected and you can ignore this.\n"
        ) + _BUILD_ADVICE + _WIKI_ADVICE
    return pct, message

def check_covariance_overlap(model_snps, covariance_snps, min_pct=1.0, sample_n=5):
    """Compares the model's SNP set against the covariance matrices' SNP set.

    Same failure shape as check_snp_overlap, one layer over: when these two
    disagree every gene is computed from zero SNPs and the results file is all
    NAs, with nothing in the log to say why. Notably, PredictDB's MASHR
    covariances are keyed by varID, so loading a MASHR model by rsid matches
    the GWAS fine and then silently fails here.
    """
    total = len(model_snps)
    if total == 0 or covariance_snps is None:
        return None, None

    overlap = model_snps & covariance_snps
    pct = 100.0 * len(overlap) / total

    model_format = SnpOverlapDiagnostics.dominant_id_format(model_snps)
    covariance_format = SnpOverlapDiagnostics.dominant_id_format(covariance_snps)
    format_mismatch = SnpOverlapDiagnostics.formats_disagree(model_format, covariance_format)

    if pct >= min_pct and not format_mismatch:
        return pct, None

    message = (
        "Only {:.2f}% of the model's SNPs were found in the covariance data ({} of {}); "
        "the model's ids look like {} and the covariance's look like {}. Every gene will "
        "be computed from zero SNPs and the results will be all NA.\n"
        "Example model SNP ids: {}\n"
        "Example covariance SNP ids: {}\n"
        "Things to check:\n"
        "- The covariance must be keyed the same way as the model: PredictDB's "
        "MASHR covariances are keyed by varID, so those models must be loaded "
        "with --model_db_snp_key varID.\n"
        "- Check that the covariance file distributed with your model is the one "
        "you passed to --covariance."
    ).format(pct, len(overlap), total, model_format, covariance_format,
             SnpOverlapDiagnostics.example_ids(model_snps, sample_n, model_format),
             SnpOverlapDiagnostics.example_ids(covariance_snps, sample_n, covariance_format))
    return pct, message

def _sanitized_gwas(gwas):
    gwas = gwas[[Constants.SNP, Constants.ZSCORE, Constants.BETA]]
    if numpy.any(~ numpy.isfinite(gwas[Constants.ZSCORE])):
        logging.warning("Discarding non finite GWAS zscores")
        gwas = gwas.loc[numpy.isfinite(gwas[Constants.ZSCORE])]
    return gwas

def _prepare_gwas(gwas):
    #If zscore is numeric, then everything is fine with us.
    # if not, try to remove "NA" strings.
    try:
        i = gwas.zscore.apply(lambda x: x != "NA")
        gwas = gwas.loc[i]
        gwas = pandas.DataFrame(gwas)
        gwas = gwas.assign(**{Constants.ZSCORE:gwas.zscore.astype(numpy.float64)})
    except Exception as e:
        logging.info("Unexpected issue preparing gwas... %s", str(e))
        pass

    if not Constants.BETA in gwas:
        gwas = gwas.assign(**{Constants.BETA: numpy.nan})

    return gwas

def _prepare_gwas_data(gwas):
    data = {}
    for x in gwas.values:
        data[x[0]] = x
    return data

def _prepare_model(model):
    K = WDBQF.K_GENE
    g = model.weights[K]
    model.weights[K] = pandas.Categorical(g, g.drop_duplicates())
    return model

def _prepare_weight_data(model, MAX_R=None):
    d,_d = [],{}
    snps = set()
    for x in model.weights.values:
        if MAX_R and len(d) + 1 > MAX_R:
            logging.info("Restricting data load to first %d", MAX_R)
            break
        gene = x[WDBQF.GENE]
        if not gene in d:
            _d[gene] = []
            d.append(gene)
        entries = _d[gene]
        entries.append(x)
        snps.add(x[WDBQF.RSID])
    return d, _d, snps

def _beta_loader(args):
    beta_contents = Utilities.contentsWithPatternsFromFolder(args.beta_folder, [])
    r = pandas.DataFrame()
    for beta_name in beta_contents:
        logging.info("Processing %s", beta_name)
        beta_path = os.path.join(args.beta_folder, beta_name)
        b = pandas.read_table(beta_path)
        r = pandas.concat([r, b])
    return r

def _gwas_wrapper(gwas):
    logging.info("Processing loaded gwas")
    return gwas

def build_context(args, gwas):
    logging.info("Loading model from: %s", args.model_db_path)
    model = load_model(args.model_db_path, args.model_db_snp_key)

    if not args.single_snp_model:
        if not args.stream_covariance:
            logging.info("Loading covariance data from: %s", args.covariance)
            covariance_manager = MatrixManager.load_matrix_manager(args.covariance)
        else:
            logging.info("Using streamed covariance from: %s", args.covariance)
            logging.warning("This version is more lenient with input covariances, as many potential errors can't be checked for the whole input covariance in advance. Pay extra care to your covariances!")
            streamer = DataFrameStreamer.data_frame_streamer(args.covariance, "GENE")
            covariance_manager = MatrixManager.StreamedMatrixManager(streamer, MatrixManager.GENE_SNP_COVARIANCE_DEFINITION)
    else:
        logging.info("Bypassing covariance for single-snp-models")
        d = model.weights[[WDBQF.K_GENE, WDBQF.K_RSID]].rename(columns={WDBQF.K_RSID:"id1"})
        d = d.assign(id2=d.id1, value=1)
        covariance_manager = MatrixManager.MatrixManager(d, {MatrixManager.K_MODEL:WDBQF.K_GENE, MatrixManager.K_ID1:"id1", MatrixManager.K_ID2:"id2", MatrixManager.K_VALUE:"value"})

    gwas = _gwas_wrapper(gwas) if gwas is not None else _beta_loader(args)
    context = _build_context(model, covariance_manager, gwas, args.MAX_R)
    return context

def _build_context(model, covariance_manager, gwas, MAX_R=None):
    gwas = _prepare_gwas(gwas)
    gwas = _sanitized_gwas(gwas)
    context = OptimizedContext(gwas, model, covariance_manager, MAX_R)
    return context

def _build_simple_context(model, covariance_manager, gwas):
    model = _prepare_model(model)
    gwas = _prepare_gwas(gwas)
    gwas = _sanitized_gwas(gwas)
    context = SimpleContext(gwas, model, covariance_manager)
    return context

def _to_int(d):
    r = d
    try:
        r = int(d)
    except:
        pass
    return r

def _results_column_order(with_additional=False):
    K = Constants
    AK = AssociationCalculation.ARF
    column_order = [WDBQF.K_GENE,
                    WDBEQF.K_GENE_NAME,
                    K.ZSCORE,
                    AK.K_EFFECT_SIZE,
                    Constants.PVALUE,
                    AK.K_VAR_G,
                    WDBEQF.K_PRED_PERF_R2,
                    WDBEQF.K_PRED_PERF_PVAL,
                    WDBEQF.K_PRED_PERF_QVAL,
                    AK.K_N_SNPS_USED,
                    AK.K_N_SNPS_IN_COV,
                    WDBEQF.K_N_SNP_IN_MODEL]
    if with_additional:
        ADD = AssociationCalculation.ASF
        column_order.extend([ADD.K_BEST_GWAS_P, ADD.K_LARGEST_WEIGHT])

    return column_order

def format_output(results, context, remove_ens_version):
    results = results.drop("n_snps_in_model",axis=1)

    # Dodge the use of cdf on non finite values
    i = numpy.isfinite(results.zscore)
    results[Constants.PVALUE] = numpy.nan
    results.loc[i, Constants.PVALUE] = 2 * stats.norm.sf(numpy.abs(results.loc[i, Constants.ZSCORE].values))

    model_info = pandas.DataFrame(context.get_model_info())

    merged = pandas.merge(results, model_info, how="inner", on="gene")
    if remove_ens_version:
        merged.gene = merged.gene.str.split(".").str.get(0)


    column_order = _results_column_order()
    merged = merged[column_order]

    # since we allow NA in covs, we massage it a little bit into resemblying an int instead of a float
    # (pandas uses the NaN float.)

    merged.n_snps_in_cov = merged.n_snps_in_cov.apply(_to_int)
    merged = merged.sort_values(by=Constants.PVALUE)
    merged = merged.fillna("NA")
    return merged

def merge_additional_output(results, stats, context, remove_ens_version):
    merged = pandas.merge(results, stats, how="inner", on="gene")
    if remove_ens_version:
        merged.gene = merged.gene.str.split(".").str.get(0)

    column_order = _results_column_order(with_additional=True)
    merged = merged[column_order]
    return merged

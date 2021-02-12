#! /usr/bin/env python
import logging
import os
from timeit import default_timer as timer

import numpy
import pandas
import gzip
import pyliftover

import metax
from metax import Utilities
from metax import Logging
from metax import Exceptions
from metax import PredictionModel
from metax.genotype import Genotype
from metax.misc import GWASAndModels, Genomics, KeyedDataSource

GF = Genotype.GF

def dosage_generator(args, variant_mapping=None, weights=None):
    if args.liftover:
        logging.info("Acquiring liftover conversion")
        liftover_chain = pyliftover.LiftOver(args.liftover)
        liftover_conversion = lambda chr,pos: Genomics.lift(liftover_chain, chr, pos, args.zero_based_positions)
    else:
        liftover_chain = None
        liftover_conversion = None

    whitelist = None
    if variant_mapping and type(variant_mapping) == dict:
        logging.info("Setting whitelist from mapping keys")
        whitelist = set(variant_mapping.keys())
    else:
        logging.info("Setting whitelist from available models")
        whitelist = set(weights.rsid)

    d = None
    if args.text_genotypes:
        from metax.genotype import DosageGenotype
        d = DosageGenotype.dosage_files_geno_lines(args.text_genotypes, variant_mapping=variant_mapping,
                whitelist=whitelist, skip_palindromic=args.skip_palindromic, liftover_conversion=liftover_conversion)
    elif args.bgen_genotypes:
        from metax.genotype import BGENGenotype
        d = BGENGenotype.bgen_files_geno_lines(args.bgen_genotypes,
            variant_mapping=variant_mapping, force_colon=args.force_colon, use_rsid=args.bgen_use_rsid, whitelist=whitelist, skip_palindromic=args.skip_palindromic)
    elif args.vcf_genotypes:
        from metax.genotype import CYVCF2Genotype
        d = CYVCF2Genotype.vcf_files_geno_lines(args.vcf_genotypes, mode=args.vcf_mode, variant_mapping=variant_mapping, whitelist=whitelist, skip_palindromic=args.skip_palindromic, liftover_conversion=liftover_conversion)

    if d is None:
        raise Exceptions.InvalidArguments("unsupported genotype input")
    if args.force_mapped_metadata:
        d = Genotype.force_mapped_metadata(d, args.force_mapped_metadata)
    return d

def model_structure(args):
    model = PredictionModel.load_model(args.model_db_path, args.model_db_snp_key)
    m = {}
    weights, extra = model.weights, model.extra
    if args.sub_batches is not None and args.sub_batch is not None:
        logging.info("slicing models")
        extra = Utilities.sub_batch(extra, args.sub_batches, args.sub_batch)
        weights = weights[weights.gene.isin(extra.gene)].reset_index(drop = True)

    if args.only_entries:
        extra = extra[extra.gene.isin(set(args.only_entries))]
        weights = weights[weights.gene.isin(set(args.only_entries))]

    for i in weights.itertuples():
        if not i.rsid in m:
            m[i.rsid] = (i.non_effect_allele, i.effect_allele, {})
        m[i.rsid][2][i.gene] = i.weight
    return m, weights, extra


def load_samples(args):
    s = None
    if args.text_sample_ids:
        if len(args.text_sample_ids) == 1:
            s = pandas.read_table(args.text_sample_ids[0], header=None, names=["FID", "IID"])
        elif args.text_sample_ids[1] == "UKB":
            k = pandas.read_table(args.text_sample_ids[0], sep=" ")
            k = k[k.sex != "D"].reset_index(drop=True)
            s = k[["ID_1", "ID_2"]].rename(columns={"ID_1": "FID", "ID_2": "IID"})
    elif args.vcf_genotypes:
        from metax.genotype import CYVCF2Genotype
        s = CYVCF2Genotype.get_samples(args.vcf_genotypes[0])
    elif args.bgen_genotypes:
        from metax.genotype import BGENGenotype
        s = BGENGenotype.get_samples(args.bgen_genotypes[0])
    elif args.generate_sample_ids:
        s = ["ID_{}".format(x) for x in range(0, args.generate_sample_ids)]
        s = [(x, x) for x in s]
        s = pandas.DataFrame(data=s, columns=["FID", "IID"])

    if s is None:
        raise Exceptions.InvalidArguments("Unsupported samples argument")
    return s

def get_variant_mapping(args, weights):
    mapping = None

    if len(args.variant_mapping):
        if len(args.variant_mapping) == 3:
            logging.info("Acquiring variant mapping")
            mapping = KeyedDataSource.load_data(args.variant_mapping[0], args.variant_mapping[1], args.variant_mapping[2], value_white_list=set(weights.rsid))
            # if args.variant_mapping[1] == "UKB":
            #     mapping = KeyedDataSource.load_data(args.variant_mapping[0], "variant", "panel_variant_id", value_white_list=set(weights.rsid))
            # elif args.variant_mapping[1] == "RSID":
            #     mapping = KeyedDataSource.load_data(args.variant_mapping[0], "variant", "rsid", value_white_list=set(weights.rsid))
            # elif args.variant_mapping[1] == "ID_TO_RSID":
            #     mapping = KeyedDataSource.load_data(args.variant_mapping[0], "id", "rsid", value_white_list=set(weights.rsid))
        else:
            raise Exceptions.InvalidArguments("Unsupported variant mapping argument")
    elif len(args.on_the_fly_mapping):
        checklist = set(weights.rsid)

    if len(args.on_the_fly_mapping) > 0:
        logging.info("Acquiring on-the-fly mapping")
        if args.on_the_fly_mapping[0] == "METADATA":
            if mapping:
                _mapping = mapping # Python scope subtlety, they are not blocks like swift
                mapping = lambda chromosome, position, ref_allele, alt_allele: Genomics.map_on_the_fly(_mapping, args.on_the_fly_mapping[1], chromosome, position, ref_allele, alt_allele)
            else:
                mapping = lambda chromosome, position, ref_allele, alt_allele: Genomics.coordinate_format(checklist, args.on_the_fly_mapping[1], chromosome, position, ref_allele, alt_allele)
        else:
            raise RuntimeError("Unsupported on_the_fly argument")
    return mapping

def prepare_prediction(args, extra, samples):
    logging.info("Preparing prediction")
    results = None
    if len(args.prediction_output) < 2:
        from metax.predixcan.Utilities import BasicPredictionRepository
        results = BasicPredictionRepository(samples, extra, args.prediction_output[0])
    else:
        if args.prediction_output[1] == "HDF5":
            from metax.predixcan.Utilities import HDF5PredictionRepository
            results = HDF5PredictionRepository(samples, extra, args.prediction_output[0])
        else:
            raise Exceptions.InvalidArguments("Unsupported output specification")
    return results

def run(args):
    start = timer()
    if args.prediction_output:
        if os.path.exists(args.prediction_output[0]):
            logging.info("Prediction output exists. Move or remove if you want this ran again.")
            return
        Utilities.ensure_requisite_folders(args.prediction_output[0])

    if args.prediction_summary_output:
        if os.path.exists(args.prediction_summary_output):
            logging.info("Summary output exists. Move or remove if you want this ran again.")
            return
        Utilities.ensure_requisite_folders(args.prediction_output[0])

    logging.info("Loading samples")
    samples = load_samples(args)

    logging.info("Loading model")
    model, weights, extra = model_structure(args)

    variant_mapping = get_variant_mapping(args, weights)

    logging.info("Preparing genotype dosages")
    dosage_source = dosage_generator(args, variant_mapping, weights)

    logging.info("Processing genotypes")
    dcapture = []
    reporter = Utilities.PercentReporter(logging.INFO, len(set(weights.rsid.values)))
    snps_found = set()
    with prepare_prediction(args, extra, samples) as results:

        for i,e in enumerate(dosage_source):
            if isinstance(e, RuntimeError):
                raise e
            if args.stop_at_variant and i>args.stop_at_variant:
                break
            var_id = e[GF.RSID]

            logging.log(8, "variant %i:%s", i, var_id)
            if var_id in model:
                s = model[var_id]
                ref_allele, alt_allele = e[GF.REF_ALLELE], e[GF.ALT_ALLELE]

                allele_align, strand_align = GWASAndModels.match_alleles(ref_allele, alt_allele, s[0], s[1])
                if not allele_align or not strand_align:
                    continue

                dosage = e[GF.FIRST_DOSAGE:]
                if allele_align == -1:
                    dosage = tuple(map(lambda x: 2 - x, dosage))
                dosage = numpy.array(dosage, dtype=numpy.float)

                snps_found.add(var_id)

                for gene, weight in s[2].items():
                    results.update(gene, dosage, weight)
                    if args.capture:
                        dcapture.append((gene, weight, var_id, s[0], s[1], ref_allele, alt_allele, strand_align, allele_align) + e[GF.FIRST_DOSAGE:])

                reporter.update(len(snps_found), "%d %% of models' snps used")

    reporter.update(len(snps_found), "%d %% of models' snps used", force=True)
     
    if args.capture:
        logging.info("Saving data capture")
        Utilities.ensure_requisite_folders(args.capture)
        with gzip.open(args.capture, "w") as f:
            header = "gene\tweight\tvariant_id\tref_allele\teff_allele\ta0\ta1\tstrand_align\tallele_align\t" + "\t".join(samples.IID.values) + "\n"
            f.write(header.encode())
            for c in dcapture:
                l = "\t".join(map(str, c)) + "\n"
                f.write(l.encode())

    if args.prediction_output and len(args.prediction_output) < 2:
        logging.info("Storing prediction")
        results.store_prediction()

    if args.prediction_summary_output:
        logging.info("Saving summary")
        summary = results.summary()
        Utilities.save_dataframe(summary, args.prediction_summary_output)

    end = timer()
    logging.info("Successfully predicted expression in %s seconds"%(str(end-start)))

    return results

def add_arguments(parser):
    parser.add_argument("--model_db_path",
                        help="Name of model db in data folder. "
                             "PrediXcan will filter input GWAS snps that are not present in the model")

    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")
    parser.add_argument("--liftover", help = "File with liftover chain. Liftover conversion will take place before any variant mapping or whitelisting.")
    parser.add_argument("--zero_based_positions", help="chromosome postions start at 0", action="store_true")
    parser.add_argument('--skip_palindromic', action="store_true", help="ignore palindromic variants (i.e. C/G)")
    parser.add_argument("--stop_at_variant", help="convenience to do an early exit", type=int, default=None)
    parser.add_argument('--bgen_genotypes', nargs='+', help="genotypes (bgen format) to use")
    parser.add_argument('--bgen_use_rsid', action="store_true", help="use rsid if available")
    parser.add_argument("--vcf_genotypes", nargs='+', help="genotypes (vcf format) to use")
    parser.add_argument("--vcf_mode", help="genotyped or imputed")
    parser.add_argument('--force_colon', action="store_true", help ="will convert variant ids from 'chr:pos_a0_a1' to 'chr:pos:a0:a1'")
    parser.add_argument('--force_mapped_metadata', help="will convert variant ids from 'chr:pos_a0_a1' to 'chr:pos:a0:a1'")
    parser.add_argument('--text_genotypes', nargs='+', help="genotypes to use")
    parser.add_argument('--text_sample_ids', help="path to file with individual samples' ids", nargs="+", default=[])
    parser.add_argument("--generate_sample_ids", help="Specify a number of samples, the ordinal number will be used as individual id", type=int, default=None)
    parser.add_argument("--prediction_output", help="name of file to put results in", nargs="+", default=[])
    parser.add_argument("--prediction_summary_output", help="name of file to put summary results in")
    parser.add_argument("--variant_mapping", help="Table to convert from genotype variants to model variants.", nargs="+", default=[])
    parser.add_argument("--on_the_fly_mapping", help="Option to convert input genotype metadata to a variant id; this can then be used with a variant mapping or directly match the models.", nargs="+", default=[])
    parser.add_argument("--sub_batches", help="split data in slices", type=int, default=None)
    parser.add_argument("--sub_batch", help="compute on a specific slice of data", type=int, default=None)
    parser.add_argument("--only_entries", help="Compute only these entries in the models (e.g. a whitelist of genes)", nargs="+")
    parser.add_argument("--capture")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Predict.py %s: Predict transcriptome from a model.' % (metax.__version__))
    add_arguments(parser)
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = 10, type = int)
    parser.add_argument("--throw", action="store_true", help="Throw exception on error")

    args = parser.parse_args()

    Logging.configureLogging(args.verbosity)

    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException as e:
            logging.error("Error:%s", e.msg)
            exit(1)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)

#! /usr/bin/env python
import logging
import os
from timeit import default_timer as timer

import numpy
import pandas

import metax
from metax import Utilities
from metax import Logging
from metax import Exceptions
from metax import PredictionModel
from metax.genotype import Genotype
from metax.misc import GWASAndModels
from metax.predixcan.Utilities import BasicPredictionRepository
from metax.data_management import KeyedDataSource

GF = Genotype.GF

def dosage_generator(args, variant_mapping=None, weights=None):
    whitelist = None
    if not variant_mapping and weights is not None:
        whitelist = set(weights.rsid)

    d = None
    if args.text_genotypes:
        from metax.genotype import DosageGenotype
        d = DosageGenotype.dosage_files_geno_lines(args.text_genotypes, snps=whitelist)
    elif args.bgen_genotypes:
        from metax.genotype import BGENGenotype
        d = BGENGenotype.bgen_files_geno_lines(args.bgen_genotypes, variant_mapping, args.force_colon, args.bgen_use_rsid, whitelist)
    if d is None:
        raise Exceptions.InvalidArguments("unsupported genotype input")
    if args.force_mapped_metadata:
        d = Genotype.force_mapped_metadata(d, args.force_mapped_metadata)
    return d

def model_structure(args):
    model = PredictionModel.load_model(args.model_db_path, args.model_db_snp_key)
    m = {}
    for i in model.weights.itertuples():
        if not i.rsid in m:
            m[i.rsid] = (i.non_effect_allele, i.effect_allele, {})
        m[i.rsid][2][i.gene] = i.weight
    return m, model.weights, model.extra


def load_samples(args):
    s = None
    if args.generate_sample_ids:
        s = ["ID_{}".format(x) for x in range(0, args.generate_sample_ids)]
        s = [(x, x) for x in s]
        s = pandas.DataFrame(data=s, columns=["FID", "IID"])
    elif len(args.text_sample_ids) == 1:
        s = pandas.read_table(args.text_sample_ids[0], header=None, names = ["FID", "IID"])
    else:
        if args.text_sample_ids[1] == "UKB":
            k = pandas.read_table(args.text_sample_ids[0], sep=" ")
            k = k[k.sex != "D"].reset_index(drop=True)
            s = k[["ID_1", "ID_2"]].rename(columns={ "ID_1": "FID", "ID_2": "IID" })
    if s is None:
        raise Exceptions.InvalidArguments("Unsupported samples argument")
    return s

def get_variant_mapping(args, weights):
    mapping = None
    if len(args.variant_mapping) == 2:
        logging.info("Acquiring variant mapping")
        if args.variant_mapping[1] == "UKB":
            mapping = KeyedDataSource.load_data(args.variant_mapping[0], "variant", "panel_variant_id", value_white_list=set(weights.rsid))
        elif args.variant_mapping[1] == "RSID":
            mapping = KeyedDataSource.load_data(args.variant_mapping[0], "variant", "rsid", value_white_list=set(weights.rsid))
        else:
            raise Exceptions.InvalidArguments("Unsupported variant mapping argument")
    return mapping

def run(args):
    start = timer()
    if args.prediction_output and os.path.exists(args.prediction_output[0]):
        logging.info("Prediction output exists. Move or remove if you want this ran again.")
        return

    if args.prediction_summary_output and os.path.exists(args.prediction_summary_output):
        logging.info("Summary output exists. Move or remove if you want this ran again.")
        return

    logging.info("Loading samples")
    samples = load_samples(args)

    logging.info("Loading model")
    model, weights, extra = model_structure(args)

    logging.info("Preparing prediction")
    if len(args.prediction_output) < 2:
        results = BasicPredictionRepository(samples, extra)

    variant_mapping = get_variant_mapping(args, weights)

    logging.info("Preparing genotype dosages")
    dosage_source = dosage_generator(args, variant_mapping, weights)

    logging.info("Processing genotypes")
    for i,e in enumerate(dosage_source):
        if args.stop_at_variant and i>args.stop_at_variant:
            break
        var_id = e[GF.RSID]
        logging.log(8, "variant %i:%s", i, var_id)
        if var_id in model:
            s = model[var_id]
            allele_align, strand_align = GWASAndModels.match_alleles(e[GF.REF_ALLELE], e[GF.ALT_ALLELE], s[0], s[1])
            if not allele_align or not strand_align:
                continue

            dosage = numpy.array(e[GF.FIRST_DOSAGE:])
            for k,v in s[2].items():
                weight = v * allele_align
                results.update(k, dosage, weight)

    if args.prediction_output and len(args.prediction_output) < 2:
        results.store_prediction(args.prediction_output[0])

    if args.prediction_summary_output:
        summary = results.summary()
        Utilities.save_dataframe(summary, args.prediction_summary_output)

    end = timer()
    logging.info("Successfully predicted expression in %s seconds"%(str(end-start)))

    return results

def add_arguments(parser):
    parser.add_argument("--model_db_path",
                        help="Name of model db in data folder. "
                             "If supplied, will filter input GWAS snps that are not present; this script will not produce output if any error is encountered."
                             "If not supplied, will convert the input GWAS as found, one line at a atime, until finishing or encountering an error.")

    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")
    parser.add_argument("--stop_at_variant", help="convenience to do an early exit", type=int, default=None)
    parser.add_argument('--bgen_genotypes', nargs='+', help="genotypes (bgen format) to use")
    parser.add_argument('--bgen_use_rsid', action="store_true", help="use rsid if available")
    parser.add_argument('--force_colon', action="store_true", help ="will convert variant ids from 'chr:pos_a0_a1' to 'chr:pos:a0:a1'")
    parser.add_argument('--force_mapped_metadata', help="will convert variant ids from 'chr:pos_a0_a1' to 'chr:pos:a0:a1'")
    parser.add_argument('--text_genotypes', nargs='+', help="genotypes to use")
    parser.add_argument('--text_sample_ids', help="path to file with individual samples' ids", nargs="+", default=[])
    parser.add_argument("--generate_sample_ids", help="Speicify a number of samples, the ordinal number will be used as individual id", type=int, default=None)
    parser.add_argument("--prediction_output", help="name of file to put results in", nargs="+", default=[])
    parser.add_argument("--prediction_summary_output", help="name of file to put summary results in")
    parser.add_argument("--variant_mapping", help="Table to convert from genotype variants to model variants", nargs="+", default=[])

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
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)
import logging
import os
from timeit import default_timer as timer

import numpy
import pandas

import metax
from metax import Utilities, Logging, Exceptions, PredictionModel
from metax.genotype import DosageGenotype, Genotype
from metax.misc import GWASAndModels
from metax.predixcan.Utilities import BasicPredictionRepository

GF = Genotype.GF

def dosage_generator(args):
    return DosageGenotype.dosage_files_geno_lines(args.text_genotypes)

def model_structure(args):
    model = PredictionModel.load_model(args.model_db_path, args.model_db_snp_key)
    m = {}
    for i in model.weights.itertuples():
        if not i.rsid in m:
            m[i.rsid] = (i.non_effect_allele, i.effect_allele, {})
        m[i.rsid][2][i.gene] = i.weight
    return m, model.weights, model.extra


def load_samples(args):
    return pandas.read_table(args.text_sample_ids, header=None, names = ["FID", "IID"])

def run(args):
    start = timer()
    if args.prediction_output and os.path.exists(args.prediction_output[0]):
        logging.info("Prediction output exists. Move or remove if you want this ran again.")
        return

    if args.prediction_summary_output and os.path.exists(args.prediction_summary_output):
        logging.info("Summary output exists. Move or remove if you want this ran again.")
        return

    dosage_source = dosage_generator(args)

    samples = load_samples(args)

    logging.info("Loading model")
    model, weights, extra = model_structure(args)

    logging.info("Preparing prediction")
    if len(args.prediction_output) < 2:
        results = BasicPredictionRepository(samples, extra)

    logging.info("Processing genotypes")
    for e in dosage_source:
        var_id = e[GF.RSID]
        if var_id in model:
            s = model[var_id]
            allele_align, strand_align = GWASAndModels.match_alleles(e[GF.REF_ALLELE], e[GF.ALT_ALLELE], s[0], s[1])
            if not allele_align or not strand_align:
                continue

            dosage = numpy.array(e[GF.FIRST_DOSAGE:])
            for k,v in s[2].iteritems():
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

    parser.add_argument('--text_genotypes', nargs='+', help="genotypes to use")
    parser.add_argument('--text_sample_ids', help="path to file with individual samples' ids")
    parser.add_argument("--prediction_output", help="name of file to put results in", nargs="+", default=[])
    parser.add_argument("--prediction_summary_output", help="name of file to put summary results in")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Predict.py %s: Predict transcriptome from a model.' % (metax.__version__))

    add_arguments(parser)
    parser.add_argument("--verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = 10, type = int)
    parser.add_argument("--throw", action="store_true", help="Throw exception on error")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))
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
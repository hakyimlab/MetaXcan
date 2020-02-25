#! /usr/bin/env python
import logging
import os
from timeit import default_timer as timer

import numpy
import pandas
import gzip

import metax
from metax import Utilities
from metax import Logging
from metax import Exceptions
from metax import PredictionModel
from metax.genotype import Genotype
from metax.misc import GWASAndModels
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
    elif args.vcf_genotypes:
        from metax.genotype import PYVCFGenotype
        d = PYVCFGenotype.vcf_files_geno_lines(args.vcf_genotypes, args.vcf_mode, whitelist=whitelist)

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
    elif args.vcf_genotypes and not args.text_sample_ids:
        from metax.genotype import PYVCFGenotype
        s = PYVCFGenotype.get_samples(args.vcf_genotypes[0])
    elif args.generate_sample_ids:
        s = ["ID_{}".format(x) for x in range(0, args.generate_sample_ids)]
        s = [(x, x) for x in s]
        s = pandas.DataFrame(data=s, columns=["FID", "IID"])

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
    print("capture: ", args.capture)
    with prepare_prediction(args, extra, samples) as results:
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

                dosage = e[GF.FIRST_DOSAGE:]
                if allele_align == -1:
                    dosage = tuple(map(lambda x: 2 - x, dosage))
                dosage = numpy.array(dosage, dtype=numpy.float)

                for gene, weight in s[2].items():
                    results.update(gene, dosage, weight)
                    if args.capture:
                        dcapture.append((gene, weight, var_id, s[0], s[1], e[GF.REF_ALLELE], e[GF.ALT_ALLELE]) + e[GF.FIRST_DOSAGE:])


    if args.capture:
        logging.info("Saving data capture")
        with gzip.open(args.capture, "w") as f:
            header = "gene\tweight\tvariant_id\tref_allele\teff_allele\ta0\ta1\t" + "\t".join( ["ID_{}".format(x) for x in range(0, len(dcapture[0][7:]))]) + "\n"
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
                             "If supplied, will filter input GWAS snps that are not present; this script will not produce output if any error is encountered."
                             "If not supplied, will convert the input GWAS as found, one line at a atime, until finishing or encountering an error.")

    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")
    parser.add_argument("--stop_at_variant", help="convenience to do an early exit", type=int, default=None)
    parser.add_argument('--bgen_genotypes', nargs='+', help="genotypes (bgen format) to use")
    parser.add_argument('--bgen_use_rsid', action="store_true", help="use rsid if available")
    parser.add_argument("--vcf_genotypes", nargs='+', help="genotypes (vcf format) to use")
    parser.add_argument("--vcf_mode", help="genotyped or imputed")
    parser.add_argument('--force_colon', action="store_true", help ="will convert variant ids from 'chr:pos_a0_a1' to 'chr:pos:a0:a1'")
    parser.add_argument('--force_mapped_metadata', help="will convert variant ids from 'chr:pos_a0_a1' to 'chr:pos:a0:a1'")
    parser.add_argument('--text_genotypes', nargs='+', help="genotypes to use")
    parser.add_argument('--text_sample_ids', help="path to file with individual samples' ids", nargs="+", default=[])
    parser.add_argument("--generate_sample_ids", help="Speicify a number of samples, the ordinal number will be used as individual id", type=int, default=None)
    parser.add_argument("--prediction_output", help="name of file to put results in", nargs="+", default=[])
    parser.add_argument("--prediction_summary_output", help="name of file to put summary results in")
    parser.add_argument("--variant_mapping", help="Table to convert from genotype variants to model variants", nargs="+", default=[])
    parser.add_argument("--sub_batches", help="split data in slices", type=int, default=None)
    parser.add_argument("--sub_batch", help="compute on a specific slice of data", type=int, default=None)
    parser.add_argument("--only_entries", help="Compute only these genes", nargs="+")
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
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)
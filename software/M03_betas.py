#! /usr/bin/env python
"""
This module contains high level code to parse GWAS/GWAMA summary statistics files and parse them into a standard format.
You can use it as stand alone tool to align it to Predictdb Models or jus tconvert the format,
or it can be called from another script to load the data into memory

TODO:
    Implement "streaming read" where the whole file is not read at once, and output as data is read from file.
    (mostly meant to make the unfiltered operation of this script have less memory footprint)

"""
__author__ = 'heroico'
import logging
import os
import re
import pandas

from timeit import default_timer as timer

import metax
__version__ = metax.__version__


from metax import Constants
from metax.misc import GWASAndModels
from metax.gwas import GWAS
from metax.gwas import Utilities as GWASUtilities
from metax import PredictionModel
from metax import Utilities
from metax import Logging
from metax import Exceptions

def build_betas(args, model, gwas_format, name, model_snp_map):
    logging.info("Building beta for %s and %s", name, args.model_db_path if args.model_db_path else "no database")

    load_from = os.path.join(args.gwas_folder, name) if args.gwas_folder else name

    snps = model.snps() if model else None
    b = GWAS.load_gwas(load_from, gwas_format, snps=snps, separator=args.separator,
            skip_until_header=args.skip_until_header, handle_empty_columns=args.handle_empty_columns, input_pvalue_fix=args.input_pvalue_fix, keep_non_rsid=args.keep_non_rsid)

    if model_snp_map:
        logging.info("Loading mapping")
        PF = PredictionModel.WDBQF
        snp_map = pandas.read_table(model_snp_map)
        snp_map_ = snp_map.rename(columns={"a0":PF.K_NON_EFFECT_ALLELE, "a1":PF.K_EFFECT_ALLELE})[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE, "panel_variant_id", "panel_variant_a0", "panel_variant_a1", "swap"]].drop_duplicates()

        logging.info("Mapping variants")
        columns = [x for x in b.columns.values]
        b = GWASAndModels.align_data_to_alleles(b, snp_map_, Constants.SNP, PF.K_RSID)
        if GWAS.ZSCORE in b:
            b = b.assign(zscore = b.zscore * b.swap)
        if GWAS.BETA in b:
            b = b.assign(beta = b.beta * b.swap)
        b = b.rename(columns={GWAS.SNP:"gwas_snp", GWAS.EFFECT_ALLELE:"gwas_effect_allele", GWAS.NON_EFFECT_ALLELE:"gwas_non_effect_allele"})\
                .drop(columns=[GWASAndModels.EA_BASE, GWASAndModels.NEA_BASE])\
                .rename(columns={"panel_variant_id":GWAS.SNP, "panel_variant_a0":GWASAndModels.NEA, "panel_variant_a1":GWASAndModels.EA})\
                [["gwas_snp", "gwas_effect_allele", "gwas_non_effect_allele"]+columns]

    if model is not None:
        logging.info("Aligning GWAS to models")
        PF = PredictionModel.WDBQF
        base = model.weights[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE]].drop_duplicates()
        b = GWASAndModels.align_data_to_alleles(b, base, Constants.SNP, PF.K_RSID)
        b = b.drop(columns=[GWASAndModels.EA_BASE, GWASAndModels.NEA_BASE])

    b = b.fillna("NA")

    if model is not None:
        logging.info("Trimming output")
        keep = [GWAS.SNP, GWAS.ZSCORE]
        if GWAS.BETA in b: keep.append(GWAS.BETA)
        b = b[keep]

    return b

def validate(args):
    if (args.gwas_file and args.gwas_folder) or (not args.gwas_file and  not args.gwas_folder):
        raise Exceptions.InvalidArguments("Provide either (--gwas_file) or (--gwas_folder [--gwas_file_pattern])")

def run(args):
    start = timer()
    validate(args)

    if args.output_folder and args.output:
        logging.info("Specify either --output_folder or --output, not both")
        return

    if args.gwas_folder:
        regexp = re.compile(args.gwas_file_pattern) if args.gwas_file_pattern else  None
        names = Utilities.contentsWithRegexpFromFolder(args.gwas_folder, regexp)
        names.sort() #cosmetic, because different filesystems/OS yield folders in different order

        if len(names) == 0:
            msg = "No GWAS files found on %s with pattern %s" % (args.gwas_folder, args.gwas_file_pattern,)
            raise Exceptions.ReportableException(msg)
    else:
        names = [args.gwas_file]

    gwas_format = GWASUtilities.gwas_format_from_args(args)
    GWAS.validate_format_basic(gwas_format)
    GWAS.validate_format_for_strict(gwas_format)
    model = PredictionModel.load_model(args.model_db_path, args.model_db_snp_key) if args.model_db_path else None

    if args.output_folder or args.output:
        if args.output_folder:
            if args.output_folder and not os.path.exists(args.output_folder):
                os.makedirs(args.output_folder)
        else:
            Utilities.ensure_requisite_folders(args.output)

        for i,name in enumerate(names):
            output_path = os.path.join(args.output_folder, name) if not args.output else args.output
            if args.output_folder or i==0:
                m = "w"
            else:
                m = "a"

            if os.path.exists(output_path):
                logging.info("%s already exists, delete it if you want it to be done again", output_path)
                continue

            b = build_betas(args, model, gwas_format, name, args.snp_map_file)
            c = "gzip" if ".gz" in output_path else None
            logging.info("Saving %s", output_path)
            b.to_csv(output_path, sep="\t", index=False, compression=c, mode=m)
        end = timer()
        logging.info("Successfully ran GWAS input processing in %s seconds" %(str(end - start)))
    else:
        r = []
        for name in names:
            b = build_betas(args, model, gwas_format, name, args.snp_map_file)
            r.append(b)
        r = pandas.concat(r)
        end = timer()
        logging.info("Successfully parsed input gwas in %s seconds"%(str(end-start)))

        return r

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M03_betas.py {}: Build betas from GWAS data as expected by MetaXcan.'.format(__version__))

    parser.add_argument("--model_db_path",
                        help="Name of model db in data folder. "
                             "If supplied, will filter input GWAS snps that are not present; this script will not produce output if any error is encountered."
                             "If not supplied, will convert the input GWAS as found, one line at a atime, until finishing or encountering an error.")

    parser.add_argument("--model_db_snp_key", help="Specify a key to use as snp_id")

    parser.add_argument("--gwas_file", help="Load a single GWAS file. (Alternative to providing a gwas_folder and gwas_file_pattern)")

    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study."
                        "If you provide this, you are likely to need to pass a --gwas_file_pattern argument value.")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).")

    parser.add_argument("--output_folder", help="name of folder to put results in")

    parser.add_argument("--output", help="name of file to put results in")

    GWASUtilities.add_gwas_arguments_to_parser(parser)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

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


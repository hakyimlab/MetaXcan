#! /usr/bin/env python
__author__ = 'heroico'
import metax
__version__ = metax.__version__
import logging
import gzip
import os
import re
import metax.KeyedDataSet as KeyedDataSet
import metax.WeightDBUtilities as WeightDBUtilities
import metax.GWASUtilities as GWASUtilities
import metax.MethodGuessing as MethodGuessing
import metax.Utilities as Utilities
import metax.Logging as Logging
import metax.Exceptions as Exceptions

class GetBetas(object):
    def __init__(self, args):
        self.weight_db_path = args.weight_db_path
        self.gwas_folder = args.gwas_folder
        self.output_folder = args.output_folder
        self.compressed_gwas = args.compressed_gwas
        self.args = args
        self.gwas_regexp = None
        if args.gwas_file_pattern:
            self.gwas_regexp = re.compile(args.gwas_file_pattern)

    def run(self):
        if self.args.weight_db_path:
            logging.info("Loading weight model")
            weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.weight_db_path)
        else:
            weight_db_logic = None

        names = Utilities.contentsWithRegexpFromFolder(self.gwas_folder, self.gwas_regexp)

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        if len(names) == 0:
            raise Exceptions.ReportableException("No GWAS files found on %s with pattern %s" %(self.gwas_folder, self.gwas_regexp.pattern,))

        for name in names:
            try:
                self.buildBetas(weight_db_logic,name)
            # This just means that there is some extra stuff inside that directory,
            # so I'm thinking we want to ignore it.
            except Exceptions.BadFilename as e:
                logging.info("Wrong file name: %s, skipping", e.msg)
                pass

    def buildBetas(self, weight_db_logic, name):
        output_path = os.path.join(self.output_folder, name)
        if not ".gz" in output_path:
            output_path += ".gz"
        if os.path.exists(output_path):
            logging.info("%s already exists, delete it if you want it to be done again", output_path)
            return

        logging.info("Building beta for %s and %s", name, self.weight_db_path if self.weight_db_path else "no database")
        input_path = os.path.join(self.gwas_folder, name)
        file_format = GWASUtilities.GWASFileFormat.fileFormatFromArgs(input_path, self.args)

        scheme = MethodGuessing.chooseGWASProcessingScheme(self.args, input_path)
        callback = MethodGuessing.chooseGWASCallback(file_format, scheme, weight_db_logic)
        if not weight_db_logic:
            GWASUtilities.loadGWASAndStream(input_path, output_path, self.compressed_gwas, self.args.separator, self.args.skip_until_header, callback)
        else:
            dosage_loader = GWASUtilities.GWASDosageFileLoader(input_path, self.compressed_gwas, self.args.separator, self.args.skip_until_header, callback)
            results, column_order = dosage_loader.load()

            # The following check is sort of redundant, as it exists in "saveSetsToCompressedFile".
            # It exists merely to provide different logging
            if len(results):
                def do_output(file, results, column_order):
                    file.write("\t".join(column_order)+"\n")
                    first = results[column_order[0]]
                    n = len(first)
                    for i in xrange(0,n):
                        line_comps = [str(results[c][i]) for c in column_order]
                        line = "%s\n" % "\t".join(line_comps)
                        file.write(line)

                with gzip.open(output_path, "wb") as file:
                    do_output(file, results, column_order)
            else:
                logging.info("No snps from the tissue model found in the GWAS file")
        logging.info("Successfully ran GWAS input processing")

def run(args):
    "Wrapper for common behavior for execution. "

    work = GetBetas(args)
    work.run()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M03_betas.py %s: Build betas from GWAS data.' % (__version__))

    parser.add_argument("--weight_db_path",
                        help="Name of weight db in data folder. "
                             "If supplied, will filter input GWAS snps that are not present; this script will not produce output if any error is encountered."
                             "If not supplied, will convert the input GWASas found, one line at a atime, until finishing or encountering an error.",
                        default=None)

    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    parser.add_argument("--output_folder",
                        help="name of folder to put results in",
                        default="intermediate/beta")

    parser.add_argument("--scheme",
                        help="Type of beta data preprocessing, optional. Options are: "
                        "'beta' (provide (beta or OR));"
                        "'beta_se' (provide (beta or OR) and standard error); "
                        "'beta_se_to_z' (provide (beta or OR) and standard error), and Zscore of beta will be output;"
                        "'z' (provide zscore of beta),"
                        " 'beta_sign_p' (sign of beta, and pvalue); beta_p (beta and pvalue)",
                        default=None)

    parser.add_argument("--or_column",
                    help="Name of column containing Odd Ratios in input files. Either 'OR_column' or 'beta_column' must be provided",
                    default=None)

    parser.add_argument("--pvalue_column",
                    help="Name of column containing p-value in input files.",
                    default=None)

    parser.add_argument("--beta_sign_column",
                    help="Name of column containing sign of beta in input files.",
                    default=None)

    parser.add_argument("--beta_column",
                    help="Name of column containing betas in input files. Either 'OR_column' or 'beta_column' must be provided",
                    default=None)

    parser.add_argument("--se_column",
                    help="Name of column containing standard error in input file.",
                    default=None)

    parser.add_argument("--beta_zscore_column",
                    help="Name of column containing beta's zscore in input file.",
                    default=None)

    parser.add_argument("--frequency_column",
                    help="Name of column containing frequency in input file",
                    default=None)

    parser.add_argument("--non_effect_allele_column",
                    help="Name of column containing non-effect allele in input file ('reference allele', if following PrediXcan format, and plink --dosage format philosophy)",
                    default="A2")

    parser.add_argument("--effect_allele_column",
                    help="Name of column containing effect (or dosage) allele in input file (dosage/effect allele)",
                    default="A1")

    parser.add_argument("--snp_column",
                    help="Name of column containing snp in input file",
                    default="SNP")

    parser.add_argument("--compressed_gwas",
                    help="Wether input files are gzip compressed files",
                    action="store_true",
                    default=False)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)

    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.",
                        default=None)

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
            logging.error(e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)


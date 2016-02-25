#!/usr/bin/env python
__author__ = 'heroico'
import os
import re
import logging
import metax.Person as Person
import metax.ThousandGenomesUtilities as ThousandGenomesUtilities
import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.Utilities as Utilities
import metax.DataSet as DataSet
import metax.Logging as Logging
import metax.Formats as Formats


class ProcessPrerequisites(object):
    """ Build data needed by posterior steps. Shouldn't be needed most than once per pipeline"""
    def __init__(self, args):
        self.dosage_folder = args.dosage_folder
        self.snp_list = args.snp_list
        self.output_folder = args.output_folder
        self.input_format = args.input_format
        self.output_format = args.output_format
        self.population_filters = args.population_filters
        self.individual_filters = [re.compile(x) for x in args.individual_filters]

        self.chromosome_in_name_regex = re.compile(args.file_pattern)

        self.samples_input = Utilities.samplesInputPath(self.dosage_folder)
        samples_name = os.path.split(self.samples_input)[1]
        self.samples_output = os.path.join(self.output_folder, samples_name)

    def run(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        self.buildPeople()

        self.buildDosages()

    def buildPeople(self):
        if os.path.exists(self.samples_output):
            logging.info("%s already exists, delete it if you want it figured out again", self.samples_output)
        else:
            if self.input_format == Formats.IMPUTE:
                Person.Person.filterSamples(self.samples_input, self.samples_output, self.population_filters, self.individual_filters)
            elif self.input_format == Formats.PrediXcan:
                Person.Person.filterSamples(self.samples_input, self.samples_output,
                                            population_filters=self.population_filters,
                                            individual_filters=self.individual_filters,
                                            row_delimiter="\t", skip_header=False)
            else:
                logging.info("Invalid input format: %s", self.input_format)
                assert False


    def buildDosages(self):
        if self.input_format == Formats.IMPUTE:
            self.processIMPUTEFiles()
        elif self.input_format == Formats.PrediXcan:
            self.processPrediXcanFiles()
        else:
            logging.info("Invalid input format: %s", self.input_format)
            assert False

    def processIMPUTEFiles(self):
        logging.info("Loading people")
        names = Utilities.hapNamesFromFolder(self.dosage_folder)
        all_people = Person.Person.allPeople(self.samples_input)
        selected_people_by_id = Person.Person.peopleByIdFromFile(self.samples_output)

        logging.info("Loading snps")
        snp_data_set = DataSet.DataSetFileUtilities.loadFromCompressedFile(self.snp_list)
        snp_dict = {rsid:True for rsid in snp_data_set.data}

        for name in names:
            output = os.path.join(self.output_folder, name)
            filter = ThousandGenomesUtilities.IMPUTEFilteredDosageFileBuilder()
            filter.base_path = self.dosage_folder
            filter.name = name
            filter.output_pattern = output
            filter.snp_dict = snp_dict
            filter.all_people = all_people
            filter.selected_people_by_id = selected_people_by_id

            if self.output_format == Formats.IMPUTE:
                filter.buildIMPUTE()
            elif self.output_format == Formats.PrediXcan:
                search = self.chromosome_in_name_regex.search(name)
                chr = search.group(1)
                filter.chromosome_name = chr
                filter.buildPrediXcan()
            else:
                logging.info("Invalid output format: %s", self.output_format)
                assert False
            #TODO: remove once it runs properly
            #break

    def processPrediXcanFiles(self):
        logging.info("Loading people")
        all_people = Person.Person.allPeople(self.samples_input, '\t', False)
        selected_people_by_id = Person.Person.peopleByIdFromFile(self.samples_output)
        logging.info("%d total people, %d selected", len(all_people), len(selected_people_by_id))

        logging.info("Loading snps")
        snp_data_set = DataSet.DataSetFileUtilities.loadFromCompressedFile(self.snp_list)
        snp_dict = {k:True for k in snp_data_set.data}
        print len(snp_dict.keys())

        contents = Utilities.contentsWithPatternsFromFolder(self.dosage_folder, ["dosage.txt.gz"])
        for content_name in contents:
            input_path = os.path.join(self.dosage_folder, content_name)
            fileBuilder = PrediXcanFormatUtilities.PrediXcanFormatFilteredFilesProcess(input_path, self.output_folder, content_name, all_people, selected_people_by_id, snp_dict)
            if self.output_format == Formats.IMPUTE:
                fileBuilder.buildIMPUTE()
            if self.output_format == Formats.PrediXcan:
                fileBuilder.buildPrediXcan()
            else:
                logging.info("Invalid output format: %s", self.output_format)
                assert False
            #TODO: remove once it runs properly
            #break

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Filters input genotype files based on specific population and SNP filters.')

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--dosage_folder",
                        help="higher level data folder",
                        default="data/1000GP_Phase3")
                        #default="data-p/dosagefiles-hapmap2")

    parser.add_argument("--snp_list",
                        help="List of valid snps",
                        default="data/hapmapSnpsCEU.list.gz")

    parser.add_argument("--output_folder",
                        help = "Folder where results will be output",
                        default="intermediate/filtered_1000GP_Phase3")
                        #default="intermediateX/filtered_dosagefiles-hapmap2")

    parser.add_argument("--file_pattern",
                        help="pattern for dosage file names, to extract chromosome part, used for (IMPUTE input format, PrediXcan output format) only",
                        default='1000GP_Phase3_(.*)')

    parser.add_argument('--population_filters', type=str, nargs='+',
                   help='Strings to filter people by."None" for no filter. Defaulted to "EUR"',
                   default=["EUR"])

    parser.add_argument('--individual_filters', type=str, nargs='+',
                   help='Patterns to filter people by. No filter by default.',
                   default=[])

    parser.add_argument('--input_format',
                   help='Input dosage files format. Valid options are: IMPUTE, PrediXcan',
                   default=Formats.IMPUTE)

    parser.add_argument('--output_format',
                   help='Output dosage files format. Valid options are: IMPUTE, PrediXcan',
                   default=Formats.PrediXcan)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    prerequisites = ProcessPrerequisites(args)
    prerequisites.run()
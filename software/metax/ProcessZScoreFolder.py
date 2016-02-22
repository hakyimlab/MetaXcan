#! /usr/bin/env python
__author__ = 'heroico'

import os
import re
import logging
import rpy2.robjects as robjects
from subprocess import  call
import Logging
import Gene

class ProcessFiles(object):
    def __init__(self, args):
        self.gene_digest_path = args.gene_digest_file
        self.input_folder = args.input_folder
        self.output_folder = args.output_folder

        self.file_regexp = None
        if args.pattern:
            self.file_regexp = re.compile(args.pattern)

    def run(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        gene_digest, gene_digest_by_name = Gene.Gene.loadFromDigest(self.gene_digest_path)

        robjects.r.source("Plots.R")
        contents = os.listdir(self.input_folder)
        for content in contents:
            if self.file_regexp and not self.file_regexp.match(content):
                logging.log(9, "Ignoring file %s", content)
                continue
            self.processContent(content, gene_digest)

        self.cleanup()

    def processContent(self, content, gene_digest):
        logging.log(9, "Processing %s", content)
        self.cookIntermediate(content, gene_digest)
        path = self.intermediatePath(content)
        data = robjects.r.export_data(path)

        comps = content.split(".csv")[0]
        comps = comps.split("TW_")

        qqunif_name = comps[0]+"__"+comps[1]+"__qqunif.png"
        qqunif_path = os.path.join(self.output_folder, qqunif_name)
        robjects.r.zscore_qqunif_from_data(data, qqunif_path)

        manhattan_name = comps[0]+"__"+comps[1]
        manhattan_path = os.path.join(self.output_folder, manhattan_name)
        robjects.r.do_manhattan_plot_from_data(data, manhattan_path)

    def cookIntermediate(self, content, gene_digest):
        path = self.inputPath(content)
        intermediate_path = self.intermediatePath(content)
        with open(path, "r") as input_file:
            with open(intermediate_path, "w") as intermediate_file:
                for i,line in enumerate(input_file):
                    if i==0:
                        header = self.buildHeader()
                        intermediate_file.write(header)
                        continue
                    comps = line.split(",")
                    if len(comps) == 4:
                        comps.insert(3,"NA")
                    gene_name = comps[0]
                    rest = ",".join(comps[1:])
                    gene = gene_digest[gene_name] if gene_name in gene_digest else None
                    base_position = gene.base_position if gene else "NA"
                    chr = gene.chromosome_name if gene else "NA"
                    output = ",".join([gene_name,chr,base_position,rest])
                    intermediate_file.write(output)

    def inputPath(self, content):
        path = os.path.join(self.input_folder, content)
        return  path

    def intermediatePath(self, content):
        path = os.path.join(self.output_folder, "intermediate")
        if not os.path.exists(path):
            os.makedirs(path)

        path = os.path.join(path, content)
        return path

    def buildHeader(self):
        header = "gene,chr,base_position,zscore,n,model_n,VAR_g\n"
        return header

    def cleanup(self):
        intermediate = self.intermediatePath("")
        command = "rm -rf "+intermediate
        call(command.split())

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build betas from GWAS data.')

    parser.add_argument("--input_folder",
                        help="name of folder with zscore files",
                        default="results_ew/results")

    parser.add_argument("--pattern",
                        help="name pattern to select files",
                        default="AMD.*")

    parser.add_argument("--output_folder",
                        help="where to output stuff",
                        default="results_ew/images_AMD")

    parser.add_argument("--gene_digest_file",
                        help="path of gene information digest",
                        default="data/gencode.v18.genes.patched_contigs.summary.protein")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = ProcessFiles(args)
    if args.throw:
        work.run()
    else:
        try:
            work.run()
        except Exception as e:
            logging.info("Unexpected error: %s", str(e))
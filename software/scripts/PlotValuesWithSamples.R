#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='Plot and compare against pure predixcan')
parser$add_argument('--zscore_file',
                    help='Path of zscore file, without extension',
                    default='results/zscores')
parser$add_argument('--predixcan_file',
                    help='path of predixcan results, without extension',
                    default='results/T1DWB_predixcan_zscore')

arguments <- parser$parse_args(commandArgs(TRUE))

# The following allows R to find sister scripts inside the same directory as
# the current script (when run as an executable program via Rscript)
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste(script.basename, "Plots.R", sep="/"))
process_zscore_file(arguments$zscore_file)
process_zscore_predixcan_file(arguments$predixcan_file)

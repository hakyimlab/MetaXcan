library(argparse)

parser <- ArgumentParser(description='Plot and compare against pure predixcan')
parser$add_argument('--zscore_file',
                    help='Path of zscore file, without extension',
                    default='results/zscores')
parser$add_argument('--predixcan_file',
                    help='path of predixcan results, without extension',
                    default='results/T1DWB_predixcan_zscore')

arguments <- parser$parse_args(commandArgs(TRUE))

source("Plots.R")
process_zscore_file(arguments$zscore_file)
process_zscore_predixcan_file(arguments$predixcan_file)

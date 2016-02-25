#! /usr/bin/env python

import os
from subprocess import call
import logging
import metax.Logging as Logging

__author__ = 'heroico'

EUR = "EUR"
EAS = "EAS"
AFR = "AFR"
SEL = "SEL"

STUDY_POP = [ EUR, EAS, AFR]
REF_POP = [ EUR, EAS, AFR]


# STUDY_POP = [ EUR, EAS, AFR, SEL]
# REF_POP = [ EUR, EAS, AFR, SEL]

# BETA = {
#     EUR: "intermediate/BETA_SIM_TGF_EUR",
#     EAS: "intermediate/BETA_SIM_TGF_EAS",
#     AFR: "intermediate/BETA_SIM_TGF_AFR",
# }
#
# COV = {
#     EUR: "intermediate/COV_DGNWB_TGF_EUR",
#     EAS: "intermediate/COV_DGNWB_TGF_EAS",
#     AFR: "intermediate/COV_DGNWB_TGF_AFR",
# }
#
# PREDIXCAN = {
#     EUR: "data/predixcan_s_dgnwb_eur.txt",
#     EAS: "data/predixcan_s_dgnwb_eas.txt",
#     AFR: "data/predixcan_s_dgnwb_afr.txt",
# }

BETA = {
    EUR: "intermediate/BETA_IGROWTH_TGF_EUR",
    EAS: "intermediate/BETA_IGROWTH_TGF_EAS",
    AFR: "intermediate/BETA_IGROWTH_TGF_AFR",
}

COV = {
    EUR: "intermediate/COV_DGNWB_TGF_EUR",
    EAS: "intermediate/COV_DGNWB_TGF_EAS",
    AFR: "intermediate/COV_DGNWB_TGF_AFR",
}

PREDIXCAN = {
    EUR: "data/predixcan_igrowth_dgnwb_eur.txt",
    EAS: "data/predixcan_igrowth_dgnwb_eas.txt",
    AFR: "data/predixcan_igrowth_dgnwb_afr.txt",
}

RESULTS = "results/results_igrowth"

#
def zscore_path(study, reference):
    output_name = "beta_%s_ref_%s_zscores.csv" % (study, reference)
    output_file = os.path.join(RESULTS, output_name)
    return output_file

def zscore(study, reference):
    command = "M04_zscores.py"
    command += " --covariance_folder " + COV[reference]
    command += " --beta_folder " + BETA[study]
    command += " --output_file " + zscore_path(study, reference)

    logging.log(9, command)
    call(command.split())

#
def predixcan_zscore_path(study, reference):
    output_name = "beta_%s_ref_%s_predixcan_zscores.csv" % (study, reference)
    output_file = os.path.join(RESULTS, output_name)
    return output_file

def post_process(study, reference):
    command = "PostProcessData.py"
    command += " --predixcan_file " + PREDIXCAN[study]
    command += " --zscore_file " + zscore_path(study, reference)
    command += " --output " + predixcan_zscore_path(study, reference)
    command += " --save_error_proxy"

    logging.log(9, command)
    call(command.split())

#
def plot(study, reference):
    command = "Rscript PlotValuesWithSamples.R"
    command += " --zscore_file " + zscore_path(study, reference).split(".csv")[0]
    command += " --predixcan_file " + predixcan_zscore_path(study, reference).split(".csv")[0]

    logging.log(9, command)
    call(command.split())

def run():
    for reference in REF_POP:
        for study in STUDY_POP:
            zscore(study, reference)
            post_process(study,reference)
            #plot(study, reference)

if __name__ == "__main__":
    Logging.configureLogging(9)
    run()

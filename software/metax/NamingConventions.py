import re
import logging
import os

re_db = re.compile(".db$")
re_tw = re.compile("^TW_")
re_0_5 = re.compile("_0.5$")
re_g_ = re.compile("^gtex_v7_")
re_i_ = re.compile("_imputed_europeans_tw_0.5_signif.db$")
def extract_model_name(path, name_pattern=None):
    p = os.path.split(path)[1]

    if name_pattern:
        r = re.compile(name_pattern)
        p = r.search(p).group(1)
        return p

    if re_db.search(p): p = re_db.sub("", p)
    if re_0_5.search(p): p = re_0_5.sub("", p)
    if re_tw.search(p): p = re_tw.sub("", p)
    if re_g_.search(p): p = re_g_.sub("", p)
    if re_i_.search(p): p = re_i_.sub("", p)
    return p

def X_treatment(tissue_tag, token, fix):
    comps = tissue_tag.split(token)
    pheno = comps[0]
    tissue = comps[1]
    tissue_tag = fix + "_" + comps[1]
    tissue_type = fix
    return pheno, tissue, tissue_tag, tissue_type

def parse_file_name(name, tissue_addition=None, tissue_model=None):
    tissue_tag = name #initial setting, will change while parsing
    tissue_type = None

    smr =  re.compile(".smr$")
    if "_elasticNet" in tissue_tag:
        tissue_tag = tissue_tag.split("_elasticNet")[0]
        tissue_model = "Elastic Net" if tissue_model == None else tissue_model
    elif "-unscaled" in tissue_tag:
        tissue_tag = tissue_tag.split("-unscaled")[0]
    elif smr.match(tissue_tag):
        tissue_tag = smr.split(tissue_tag)[0]
    elif ".meta.txt" in tissue_tag:
        tissue_tag = tissue_tag.split(".meta.txt")[0]
    elif ".zscores.csv" in tissue_tag:
        tissue_tag = tissue_tag.split(".zscores.csv")[0]
    elif ".csv" in tissue_tag:
        tissue_tag = tissue_tag.split(".csv")[0]
    else:
        logging.info("Bad tag %s", tissue_tag)
        return


#Parse name
    def PRS_treatment(tissue_tag, token):
        return X_treatment(tissue_tag, token, "PRS")

    def TW_treatment(tissue_tag, token):
        return X_treatment(tissue_tag, token, "TW")

    def TS_treatment(tissue_tag, token):
        return X_treatment(tissue_tag, token, "TS")

    def eQTL_Treatment(tissue_tag, token):
        return X_treatment(tissue_tag, token, "eQTL")


    def DGN_treatment(tissue_tag, token):
        return token_treatment(tissue_tag, token, "DGN_WB")

    def CrossTissue_treatment(tissue_tag, token):
        return token_treatment(tissue_tag, token, "CrossTissue")

    def token_treatment(tissue_tag, token, label):
        comps = tissue_tag.split(token)
        pheno = comps[0]
        tissue = label
        tissue_tag  = label
        tissue_type = None
        return pheno, tissue, tissue_tag, tissue_type

    if "_TW_Whole_Blood_DGN" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = DGN_treatment(tissue_tag, "_TW_Whole_Blood_DGN")
    elif "_TW_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = TW_treatment(tissue_tag, "_TW_")
    elif "-TW_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = TW_treatment(tissue_tag, "-TW_")
    elif "TW_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = TW_treatment(tissue_tag, "TW_")
    elif "_TS_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = TS_treatment(tissue_tag, "_TS_")
    elif "TS_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = TS_treatment(tissue_tag, "TS_")
    elif "_DGN" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = DGN_treatment(tissue_tag, "_DGN")
    elif "DGN" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = DGN_treatment(tissue_tag, "DGN")
    elif "_CrossTissue_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = CrossTissue_treatment(tissue_tag, "CrossTissue_")
    elif "CrossTissue_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = CrossTissue_treatment(tissue_tag, "CrossTissue_")
    elif "CrossTissue" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = CrossTissue_treatment(tissue_tag, "CrossTissue")
    elif "_PRS_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = PRS_treatment(tissue_tag, "_PRS_")
    elif "_PRS" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = PRS_treatment(tissue_tag, "_PRS")
    elif "PRS_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = PRS_treatment(tissue_tag, "PRS_")
    elif "_Intron_geuvadis" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = token_treatment(tissue_tag, "_Intron_geuvadis", "Intron-Geuvadis")
    elif "_gEUVADIS" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = token_treatment(tissue_tag, "_gEUVADIS", "Geuvadis")
    elif "_microRNA" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = token_treatment(tissue_tag, "_microRNA", "microRNA")
    elif "_eQTL_" in tissue_tag:
        pheno, tissue, tissue_tag, tissue_type = eQTL_Treatment(tissue_tag, "_eQTL_")
    else:
        logging.info("Bad name: %s %s", name, tissue_tag)
        return

    if pheno.endswith(("_")):
        pheno = pheno[:-1]

    if tissue_addition:
        tissue_tag = tissue_tag + tissue_addition

    return pheno, tissue, tissue_tag, tissue_type
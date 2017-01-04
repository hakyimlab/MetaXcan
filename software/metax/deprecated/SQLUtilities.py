#! /usr/bin/env python

import os
import psycopg2
import scipy.stats as stats
import logging
import Logging

class CSVTF1(object):
    GENE=0
    GENE_NAME=1
    ZSCORE=2
    PVALUE=3
    PRED_PERF_R2=4
    VAR_G=5
    N=6
    COVARIANCE_N=7
    MODEL_N=8

    header="gene,gene_name,zscore,pvalue,pred_perf_R2,VAR_g,n,covariance_n,model_n"

class CSVTF2(object):
    GENE=0
    GENE_NAME=1
    ZSCORE=2
    VAR_G=3
    N=4
    COVARIANCE_N=5
    MODEL_N=6
    PRED_PERF_R2=7

    header="gene,gene_name,zscore,VAR_g,n,covariance_n,model_n,pred_perf_R2"

def check_table(conn, table_name):
    if not str.isalnum(table_name):
        raise RuntimeError("Won't accept %s" %(table_name,))
    cursor = conn.cursor()
    query = 'CREATE TABLE IF NOT EXISTS metaxcanresults (' \
    ' "gene" varchar,' \
    ' "gene_name" varchar,' \
    ' "zscore" real,' \
    ' "n" integer,' \
    ' "model_n" integer,' \
    ' "pred_perf_r2" real,' \
    ' "pval" double precision,' \
    ' "tissue" varchar,' \
    ' "pheno" varchar' \
    ' );'

    cursor.execute(query)
    conn.commit()


def process_results_file(path, conn, table_name, tissue_tag=None):
    if not str.isalnum(table_name):
        raise RuntimeError("Cannot accept %s" % (table_name, ))

    content = os.path.basename(os.path.normpath(path))

    tissue = content
    if "_elasticNet" in tissue: tissue =tissue.split("_elasticNet")[0]
    if "-unscaled" in tissue: tissue = tissue.split("-unscaled")[0]
    if ".csv" in tissue: tissue = tissue.split(".csv")[0]

    if "_DGN" in tissue:
        pheno = tissue.split("_DGN")[0]
        tissue = "DGN_WB"
    elif "DGN" in tissue:
        pheno = tissue.split("DGN")[0]
        tissue = "DGN_WB"
    elif "_TW_" in tissue:
        pheno = tissue.split("_TW_")[0]
        tissue = "TW_"+tissue.split("_TW_")[1]
    elif "TW_" in tissue:
        pheno = tissue.split("TW_")[0]
        tissue = "TW_"+tissue.split("TW_")[1]
    elif "_CrossTissue_" in tissue:
        pheno = tissue.split("_CrossTissue_")[0]
        tissue = "CrossTissue"
    elif "CrossTissue_" in tissue:
        pheno = tissue.split("CrossTissue_")[0]
        tissue = "CrossTissue"
    else:
        logging.info("Bad name: %s", content)
        return

    if tissue_tag:
        tissue = tissue + tissue_tag

    logging.info("opening %s",path)
    cursor = conn.cursor()
    data = []
    with open(path) as file:
        for i, line in enumerate(file):
            if i==0:
                header = line.strip()
                if header == CSVTF1.header:
                    logging.info("selected new format")
                    RTF=CSVTF1
                elif header == CSVTF2.header:
                    logging.info("selected old format")
                    RTF=CSVTF2
                else:
                    raise RuntimeError("Invalid header")
                continue
            # if i % 1000 == 0:
            #     logging.info("Rolling!")


            comps = line.strip().split(",")
            if comps[RTF.ZSCORE] == "NA": continue
            if "inf" in comps[RTF.ZSCORE]: continue

            gene = comps[RTF.GENE]
            gene_name = comps[RTF.GENE_NAME]
            zscore = float(comps[RTF.ZSCORE])
            n = int(comps[RTF.N])
            model_n = int(comps[RTF.MODEL_N])
            pred_perf_r2 = float(comps[RTF.PRED_PERF_R2])
            if RTF == CSVTF1:
                pval = float(comps[RTF.PVALUE])
            else:
                pval = stats.norm.sf(abs(zscore)) * 2

            data.append((gene, gene_name, zscore, n, model_n, pred_perf_r2, pval, tissue, pheno, ))
            # cursor.execute(
            #         "insert into " + table_name + " (gene, gene_name, zscore, n, model_n, pred_perf_r2, pval, pheno, tissue) values (%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            #         (gene, gene_name, zscore, n, model_n, pred_perf_r2, pval, pheno, tissue,))
    logging.info("Making query")
    dataText = ','.join(cursor.mogrify('(%s,%s,%s,%s,%s,%s,%s,%s,%s)', row) for row in data)
    cursor.execute('insert into ' + table_name + ' values ' + dataText)
    conn.commit()

class LoadIntoSQL(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        logging.info("Connecting...")
        conn = self.build_connection()

        logging.info("Checking table...")
        check_table(conn,self.args.table_name)

        logging.info("processing")
        contents = os.listdir(self.args.results_folder)

        for content in contents:
            skip = False
            for e in self.args.exclude_patterns:
                if e in content:
                    skip = True
                    break
            if skip:
                logging.info("skipping %s", content)
                continue

            logging.info("processing %s", content)
            path = os.path.join(self.args.results_folder, content)
            process_results_file(path, conn, self.args.table_name, self.args.tissue_tag)

    def build_connection(self):
        conn = psycopg2.connect(host=self.args.host, database=self.args.db_name, user=self.args.user_name, password=self.args.user_password)
        return conn

if __name__== "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Load results into SQL')

    parser.add_argument("--host",
                        help="db host",
                        default="127.0.0.1")

    parser.add_argument("--db_name",
                        help="database name",
                        default="kk2")

    parser.add_argument("--user_name",
                        help="database user name")

    parser.add_argument("--user_password",
                        help="database user name")

    parser.add_argument("--table_name",
                        help="alphanumeric table name",
                        default="metaxcanresults")

    parser.add_argument('--exclude_patterns', type=str, nargs='+',
                help='Strings to filter out results files. (not regexp at the moment)',
                default=["TS_", "Organ_"])

    parser.add_argument("--results_folder",
                        help="path to folder with metaxcan results",
                        default="results")

    parser.add_argument("--tissue_tag",
                        help="String addition to tissue name",
                        default=None)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = LoadIntoSQL(args)
    work.run()
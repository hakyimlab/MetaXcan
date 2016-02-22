#! /usr/bin/env python
__author__ = 'heroico'

import logging
import os
import Gene
import Logging
import KeyedDataSet

EXCLUDED = { "MROH7":True, "DCAF8":True, "F11R":True,"LIMS3":True, "CRYBG3":True, "FAM47E":True, "CKS1B":True, "DEFB130":True, "TMEM236":True,
             "LSP1":True, "CCDC177":True, "GOLGA6L9":True, "ZNF177":True, "ZNF763":True, "PLCXD1":True, "GTPBP6":True, "PPP2R3B":True,
             "SHOX":True, "CRLF2":True, "CSF2RA":True, "IL3RA":True, "SLC25A6":True, "ASMTL":True, "P2RY8":True, "AKAP17A":True, "ASMT":True,
             "DHRSX":True, "ZBED1":True, "CD99":True, "SPRY3":True, "VAMP7":True, "IL9R":True}

M_R_HEADER = "gene,gene_name,zscore,pvalue,pred_perf_R2,VAR_g,n,covariance_n,model_n"

class PostProcessData(object):
    def __init__(self, args):
        self.gene_digest_file = args.gene_digest_file
        self.predixcan_path = args.predixcan_file
        self.zscore_path = args.zscore_file
        self.exclude = args.exclude

        self.output_path = args.output

        self.save_error_proxy = args.save_error_proxy

    def run(self):
        # this will go through the same files again. Inefficient but convenient code-wise
        gene_digest, gene_digest_by_name = Gene.Gene.loadFromDigest(self.gene_digest_file)

        gene_name = KeyedDataSet.KeyedDataSetFileUtilities.loadFromFile(self.zscore_path, ",", 1, header=M_R_HEADER)
        zscore = KeyedDataSet.KeyedDataSetFileUtilities.loadFromFile(self.zscore_path, ",", 2, header=M_R_HEADER)
        n_snp = KeyedDataSet.KeyedDataSetFileUtilities.loadFromFile(self.zscore_path, ",", 7, header=M_R_HEADER)

        predixcan_header = "gene beta z-stat p-val" if not self.save_error_proxy else "gene beta z-stat p-val se(beta)"
        predixcan = KeyedDataSet.KeyedDataSetFileUtilities.loadFromFile(self.predixcan_path, " ", 2, header=predixcan_header)

        # this will go through the same files again. Inefficient but convenient code-wise
        var_g = KeyedDataSet.KeyedDataSetFileUtilities.loadFromFile(self.zscore_path, ",", 5, M_R_HEADER) if self.save_error_proxy else None
        predixcan_beta_error = KeyedDataSet.KeyedDataSetFileUtilities.loadFromFile(self.predixcan_path, " ", 4,header=predixcan_header) if self.save_error_proxy else None

        with open(self.output_path, 'w') as file:
            lines = []
            for gene, zscore_value in zscore.values_by_key.iteritems():
                if self.exclude and gene in EXCLUDED:
                    continue
                line = self.buildLine(gene, gene_name, zscore_value, n_snp, predixcan, var_g, predixcan_beta_error, gene_digest, gene_digest_by_name)
                if line:
                    lines.append(line)
                else:
                    logging.info("1p: Discarding line for %s", gene)


            for gene, predixcan_value in predixcan.values_by_key.iteritems():
                break
                #R converts the gene names because its too dumb to handle "-" as a string.
                un_r_gene = gene.replace(".", "-")
                if not un_r_gene in zscore.values_by_key:
                    if self.exclude and gene in EXCLUDED:
                        continue
                    line = self.buildLine(un_r_gene, gene_name, "NA", n_snp, predixcan, var_g, predixcan_beta_error, gene_digest, gene_digest_by_name)
                    if line:
                        lines.append(line)
                    else:
                        logging.info("2p: discarding line for %s", un_r_gene)
            lines.sort(key=lambda l: -abs(float(l.split(",")[4])) if "NA" not in l else 0)
            header = self.buildHeader()
            file.write(header)
            for line in lines:
                file.write(line)

    def buildHeader(self):
        header = "chr,base_position,gene,gene_name,zscore,predixcan_result,n_snp\n"
        if self.save_error_proxy:
            header = "chr,base_position,gene,gene_name,zscore,predixcan_result,n_snp,VAR_g,se_predixcan_beta\n"
        return header

    def buildLine(self, gene, gene_name, zscore_value, n_snp, predixcan, var_g, predixcan_beta_error, gene_digest, gene_digest_by_name):
        if gene in gene_digest:
            digest = gene_digest[gene]
        else:
            digest = gene_digest_by_name[gene] if gene in gene_digest_by_name else None

        if not digest:
            return None

        the_gene_name = gene_name.values_by_key[gene] if gene in gene_name.values_by_key else gene
        predixcan_value = "NA"
        #R converts the gene names because its too dumb to handle "-" as a string.
        r_gene = gene.replace("-", ".")
        if not r_gene in predixcan.values_by_key:
            t = the_gene_name.replace("-", ".")
            if t in predixcan.values_by_key:
                r_gene = t

        if r_gene in predixcan.values_by_key:
            predixcan_value = predixcan.values_by_key[r_gene]

        n_snp_value = n_snp.values_by_key[gene] if gene in n_snp.values_by_key else  "NA"

        v = var_g.values_by_key[gene] if (self.save_error_proxy and gene in var_g.values_by_key) else "NA"

        p = predixcan_beta_error.values_by_key[r_gene] if (self.save_error_proxy and r_gene in predixcan_beta_error.values_by_key) else  "NA"

        base_position = digest.base_position if digest else "NA"

        chr = digest.chromosome_name if digest else "NA"

        if self.save_error_proxy:
            line = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (chr, base_position, gene, the_gene_name, zscore_value, predixcan_value, n_snp_value, v, p)
        else:
            line = "%s,%s,%s,%s,%s,%s,%s\n" % (chr, base_position, gene, the_gene_name, zscore_value, predixcan_value, n_snp_value)
        return line


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Post process results.')

    parser.add_argument("--exclude",
                        help="exclude some genes from a mystical list.",
                        action="store_true",
                        default=False)

    parser.add_argument("--gene_digest_file",
                        help="path of gene information digest",
                        default="data/gencode.v18.genes.patched_contigs.summary.protein")

    parser.add_argument("--predixcan_file",
                        help="name of PrediXcan results in data folder",
                        default="data/PrediXcan_T1D_DGNWholeBlood_EN0.5.txt")

    parser.add_argument("--zscore_file",
                        help="File with zscore results in results folder",
                        default='results/zscores.csv')

    parser.add_argument("--output",
                        help="File with merged zscore and predixcan results",
                        default="results/T1DWB_predixcan_zscore.csv")

    parser.add_argument("--save_error_proxy",
                    help="Output beta error proxy",
                    action="store_true",
                    default=False)

    args = parser.parse_args()

    Logging.configureLogging(logging.INFO)

    work = PostProcessData(args)
    work.run()
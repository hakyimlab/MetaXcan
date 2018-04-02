#!/usr/bin/env python
import os
import gzip
import sqlite3
import shutil

selected=["ENSG00000107937.14", "ENSG00000107959.11", "ENSG00000234745.5"]

def copy_covariance(ip, op):
    with gzip.open(ip) as infile:
        header = infile.readline()
        with gzip.open(op, "w") as outfile:
            outfile.write(header)
            for line in infile:
                gene = line.split()[0]
                if not gene in selected:
                    continue
                outfile.write(line)

if not os.path.exists("meta_covariance"): os.makedirs("meta_covariance")
MCP="/home/numa/Documents/Projects/data/metaxcan/data/GTExCovariance/results/snp_covariance.txt.gz"
copy_covariance(MCP,"meta_covariance/snps_covariance.txt.gz")


PMP = "/home/numa/Documents/Projects/data/metaxcan/GTEx-V6p-HapMap-2016-09-08"
if not os.path.exists("dbs_3"): os.makedirs("dbs_3")
for m in os.listdir(PMP):
    if not ".db" in m: continue
    im = os.path.join(PMP, m)
    om = os.path.join("dbs_3", m)
    shutil.copyfile(im, om)

    conn = sqlite3.connect(om)
    c = conn.cursor()
    c.execute("DELETE FROM extra WHERE gene NOT IN ({})".format(",".join(["'{}'".format(x) for x in selected])))
    c.execute("DELETE FROM weights WHERE gene NOT IN ({})".format(",".join(["'{}'".format(x) for x in selected])))
    c.execute("VACUUM")
    conn.commit()
    conn.close()

if not os.path.exists("cov_3"): os.makedirs("cov_3")
for m in os.listdir(PMP):
    if not ".txt.gz" in m: continue
    ip = os.path.join(PMP, m)
    op = os.path.join("cov_3", m)
    copy_covariance(ip, op)
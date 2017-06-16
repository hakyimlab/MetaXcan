#!/usr/bin/env python
import gzip
import numpy
snps = {'1_16393357_T_G_b37', '1_15724033_A_G_b37'}

d = {}
PATH = "/run/user/1000/gvfs/smb-share:server=bulkstorage.uchicago.edu,share=im-lab/nas40t2/tempo/genotypes/Muscle_Skeletal_Analysis.snps.txt.gz"
with gzip.open(PATH) as file:
    for l in file:
        comps = l.strip().split()
        snp = comps[0]
        if not snp in snps: continue
        d[snp] = numpy.array(comps[1:], dtype=numpy.float64)
        if len(d) == len(snps): break

import numpy

k = [v for k,v in d.iteritems()]
from IPython import embed;
embed()

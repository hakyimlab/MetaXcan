#!/usr/bin/env python

import sys
import os
import io
import gzip
import metax.Gene

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest

class GeneEntry(object):
    def __init__(self, chr, start, end, ensid, name):
        self.chr = chr
        self.start = start
        self.end = end
        self.ensid = ensid
        self.name = name

    def __str__(self):
        return "\t".join([self.chr, "0", str(self.start), str(self.end), self.ensid, self.name, '0', '0'])

class TestGenes(unittest.TestCase):
    def setUp(self):
        self.txt_filename = "__gene_data.txt"
        self.header_filename = "__gene_header.txt"
        self.gz_filename = "__gene_data.txt.gz"

        self.genes = [
            GeneEntry("chr11", '87132949', '87342639', "ENSG00000085563.14", "ABCB1"),
            GeneEntry("chr6", '32812986', '32821755', "ENSG00000168394.10", "PSF1"),
            GeneEntry("chr19", '45409011', '45412650', "ENSG00000130203.9", "APO-E"),
            GeneEntry("chr6", '32840717', '32844703', "ENSG00000204264.8", "PSMB8"),
            GeneEntry("chr15", '84199311', '94230136', "ENSG00000225151.10", "GOLGA2P7"),
            GeneEntry("chrX", '148500619', '149000663', "ENSG00000155966.13", "AFF2")
        ]
        self.by_ensembl = {}
        self.by_name = {}
        for gene in self.genes:
            self.by_ensembl[gene.ensid] = gene
            self.by_name[gene.name] = gene


        with open(self.txt_filename, 'w') as file:
            for gene in self.genes:
                print(gene, file=file)


        self.header = "CHR\tSTART\tEND\tENS\tNAME\tASDF\tFDSA"
        with open(self.header_filename, 'w') as file:
            print(self.header, file=file)
            for gene in self.genes:
                print(gene, file=file)

        with io.TextIOWrapper(gzip.open(self.gz_filename, 'w'), newline="") as file:
            for gene in self.genes:
                file.write(str(gene)+"\n")

    def tearDown(self):
        os.remove(self.txt_filename)
        os.remove(self.header_filename)
        os.remove(self.gz_filename)


    def testGeneUncompressed(self):
        gene_by_ensembl, gene_by_name = metax.Gene.Gene.loadFromDigest(self.txt_filename, compressed=False, autosome_only=False)

        self.assertEqual(len(self.genes), len(gene_by_name))
        for ensid in gene_by_ensembl:
            local_gene = self.by_ensembl[ensid]
            gene = gene_by_ensembl[ensid]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)

        for name in gene_by_name:
            local_gene = self.by_name[name]
            gene = gene_by_name[name]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)

    def testAutoOnly(self):
        gene_by_ensembl, gene_by_name = metax.Gene.Gene.loadFromDigest(self.txt_filename, compressed=False, autosome_only=True)

        self.assertEqual(len(self.genes) - 1, len(gene_by_name))
        for ensid in gene_by_ensembl:
            local_gene = self.by_ensembl[ensid]
            gene = gene_by_ensembl[ensid]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)

        for name in gene_by_name:
            local_gene = self.by_name[name]
            gene = gene_by_name[name]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)

    def testHeader(self):

        gene_by_ensembl, gene_by_name = metax.Gene.Gene.loadFromDigest(self.header_filename, header=self.header, compressed=False, autosome_only=True)

        self.assertEqual(len(self.genes) - 1, len(gene_by_name))
        for ensid in gene_by_ensembl:
            local_gene = self.by_ensembl[ensid]
            gene = gene_by_ensembl[ensid]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)

        for name in gene_by_name:
            local_gene = self.by_name[name]
            gene = gene_by_name[name]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)


    def testGeneCompressed(self):
        gene_by_ensembl, gene_by_name = metax.Gene.Gene.loadFromDigest(self.gz_filename, compressed=True, autosome_only=False)

        self.assertEqual(len(self.genes), len(gene_by_name))
        for ensid in gene_by_ensembl:
            local_gene = self.by_ensembl[ensid]
            gene = gene_by_ensembl[ensid]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)

        for name in gene_by_name:
            local_gene = self.by_name[name]
            gene = gene_by_name[name]
            self.assertEqual(local_gene.chr.replace("chr", ""), gene.chromosome_name)
            self.assertEqual(local_gene.start, gene.base_position)
            self.assertEqual(local_gene.name, gene.name)
            self.assertEqual(local_gene.ensid, gene.ens_id)


if __name__ == "__main__":
    unittest.main()

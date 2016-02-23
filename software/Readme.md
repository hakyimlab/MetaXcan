# MetaXcan

MetaXcan: summary statistics based gene-level association test

# Software Description

This repository contains software for implementing and/or using MetaXcan method.
Most of it is implemented at [Python](https://www.python.org/) scripts, 
although there are some [R](https://www.r-project.org/) scripts for processing and plotting results.

Most scripts were written as reusable components in a software library,
aimed at building different processing workflows.

Some of them, however, can be executed as command line tools to perform MetaXcan computations.

Check the root [Readme](https://github.com/hakyimlab/MetaXcan) for a sample script usage.

## Prerequisites

The software is developed and tested in Linux and Max OS environments. Should be mostly working on Windows.

You need [Python 2.7](https://www.python.org/) and [numpy](http://www.numpy.org/) to run MetaXcan.
Some support scripts use [scipy](http://www.scipy.org/) too.

[R](https://www.r-project.org/) with [ggplot](http://ggplot2.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) 
is needed for some optional statistics and charts.

# Overview

MetaXcan is concerned with obtaining gene-level association tests from ordinary GWAS data.

Ordinarily, a user would need to obtain/download support data sets comprising of:
- A Transcriptome Prediction Model database (an example is [here](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/data/DGN-WB_0.5.db))
- A file with the covariance matrices of the SNPs within each gene model (such as [this one](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/intermediate/covariance.txt.gz))

And use them to run MetaXcan analysis on:
- GWAS results (such as [these](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/data/GWAS.tar.gz))

However, if you have access to interesting data,
you can build your own Transcriptome Prediction Model database 
(using tools such as [PredictDB](https://github.com/hakyimlab/PrediXmod/tree/master/PredictDB)), 
and your own covariance (using tools from this repository).

More detailed information can be found at the [Wiki](https://github.com/hakyimlab/MetaXcan/wiki),
where more specific coverage is given to the tools contained here.

# Brief Tour

Scripts **M00_prerequisites.py**, **M01_covariances_correlations.py**, **M02_variances.py**, **M03_betas.py**
and **M04_zscores.py** are steps in a MetaXcan pipeline, although tipically a user will only need to
use the last three. For ease of use, there is a script that performs these last two steps at once, **MetaXcan.py**.
There is also a GUI app, **MetaXcanUI.py**, that offers a simplified operation for **MetaXcan.py**'s functionality.

All of these runable scripts take different number of command line parameters. Run them with
**--help** option to see the options.

## M00_prerequisites.py

This script takes extensive genotype data sets such as Thousand Genomes Haplotypes 
([here](https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html))
allows for some filtering based on population or ID, and selects those entries
present in hapmap2 biallelic SNP's.

It supports input and output into a custom format that we call internally "PrediXcan Format". It is composed of
gzip compressed textfiles, without headers, where each row contains:

```bash
chromosome SNP_id SNP_position reference_allele effect_allele allele_frequency ...
# "..." stands for allele dosage data, one value per sample individual, value in [0, 2]
```

It is assumed that a "samples file" is provided, describing each sample individual's metadata, which is a plain text file
that looks like:

```bash
ID POP GROUP SEX
HG00096 GBR EUR male
HG00097 GBR EUR female
HG00099 GBR EUR female
HG00100 GBR EUR female
HG00101 GBR EUR male
HG00102 GBR EUR female
...
```
This script is hardly ever necessary. You might use it occasionally to reduce file sizes
used at the covariance script step (**M01_covariances_correlations.py**).

## M01_covariances_correlations.py

This script builds the covariance matrices needed at **M04_zscores.py**.
You will run it once in a while, if ever.
It takes input from **M00_prerequisites.py**'s output, and a genetic expression model database, such as
[this one](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/data/DGN-WB_0.5.db).

It will build the correlation matrix between SNP's in a same gene's model, for each gene, and save them
in a gzip-compressed text file.

## M03_betas.py

This is an optional script that takes GWAS results and converts them into a format that **M04_zscores.py**
can understand, and filters down snps by removing those excluded from the genetic expression.
If you use tools such as [Plink](https://www.cog-genomics.org/plink2), you will need to do this.

This supports both compressed and uncompressed text files as input. When running this script,
you need to specify the column names of expected data. We recommend files that include SNP's GWAS results for
(*beta* or *odd ratio* or *sign of beta*) and p-value. So that, for example, you provide **-or_column** and **--pvalue_column**
when running it.

## M04_zscores.py

This performs the actual MetaXcan calculation. It needs a genetic expression model database,
the covariance matrices from this model measured on a specific population, and GWAS data
as in the output of **M03_betas.py**.

If you have the model database and covariance matrix, you can just call this at the end
of your gwas pipeline (optionally using **M03_betas.py**).

Its output is a CSV file that looks like:

```
gene,gene_name,zscore,pvalue,pred_perf_R2,VAR_g,n,covariance_n,model_n
CRIM1,CRIM1,-4.19069760413,2.78098095601e-05,0.13320775358,0.0983344808163,37,37,37
FEZ2,FEZ2,3.97075348924,7.16456810029e-05,0.261426523524,0.23453565545,24,24,24
ABLIM2,ABLIM2,-3.65940334247,0.000252803173193,0.0219411680431,0.0183096836859,30,31,31
```
Where each row is a gene's association result:
* gene: a gene's id: some Tissue models provide Ensemble Id, some others (mainly DGN Whole Blood) provide [Genquant](http://www.gencodegenes.org/)'s gene name
* gene_name: gene name extracted from Genquant
* zscore: MetaXcan'as association result for the gene
* pvalue: P-value of the aforementioned statistic.
* pred_perf_R2: R2 of tissue model's correlation to gene's measured transcriptome
* n: number of snps from GWAS that got used in MetaXcan analysis
* covariance_n: number of snps in the covariance matrix
* model_n: number of snps in the model
* VAR_g: variance of the gene expression, calculated as *W' * G * W* 
(where *W* is the vector of SNP weights in a gene's model,
*W'* is its transpose, and *G* is the covariance matrix)

## MetaXcanUI.py

This script launches a desktop app that processes GWASinput and produces MetaXcan association results.
It is basically a friendlier way to use functionality at **M03_betas.py** and **M04_zscores.py**.
It doesn't support all of the options these scripts do (that as a matter of fact haven't been covered in this readme)
but might be friendlier to non technical users.

Here is a screenshot:
![screen shot](../manuscript/plots/gui.png)

## The Rest

The other scripts in this repository fall into two categories:
* components to be used by runnable scripts (such as the ones just described) 
* Support data processing and analysis

# What do I do with this stuff?

It depends on your needs and means.

Ordinarily, most users would need to download some genetic expression databases,
some covariance matrices, and run MetaXcan on GWAS studies.
The **GUI**,  **M03_betas.py**, **M04_zscores.py** and **MetaXcan.py**
can be used for that.

A more sophisticated application could include using MetaXcan
with custom expression model databases, and different population's covariance matrices.
**M00_prerequisites.py** and **M01_covariances_correlations.py** can help with that.

A yet more sophisticated use would imply using and modifying the source code
to build custom processing pipelines.

If curious, check the [Wiki](https://github.com/hakyimlab/MetaXcan/wiki) for further information.


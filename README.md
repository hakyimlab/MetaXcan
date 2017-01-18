# MetaXcan

MetaXcan: summary statistics based gene-level association test

## Introduction

MetaXcan is an extension of [PrediXcan](https://github.com/hakyimlab/PrediXcan) method, that infers the results of PrediXcan using only summary statistics.
This repository contains a wide set of tools for calculating the gene-level association test results,
and building processing pipelines as well.

## Prerequisites

The software is developed and tested in Linux and Max OS environments. The main MetaXcan script is supported in Windows.

You need [Python 2.7](https://www.python.org/), [numpy](http://www.numpy.org/), and [scipy](http://www.scipy.org/) to run MetaXcan, and there is a GUI done in TKInter.

[R](https://www.r-project.org/) with [ggplot](http://ggplot2.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) 
is needed for some optional statistics and charts.

## Project Layout

You will find a preliminary version of MetaXcan's manuscript under **manuscript** folder.

**software** folder contains an implementation of MetaXcan's method. 
The following scripts from that folder are different steps in the MetaXcan pipeline:

```bash
M00_prerequisites.py
M01_covariances_correlations.py
M02_variances.py
M03_betas.py
M04_zscores.py
```
, although a typical user will use ony the last two of them.

The rest of the scripts in **software** folder are python packaging support scripts,
and convenience wrappers such as the GUI.

Subfolder **software/metax** contains the bulk of Metaxcan's logic, implemented as a python package.


## Input data
MetaXcan will calculate the association results from GWAS results, as output by [plink](https://www.cog-genomics.org/plink2).
Some support data is needed, that needs to be set up prior MetaXcan execution.

The gist of MetaXcan input is:
- A Transcriptome Prediction Model database (an example is [here](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/DGN-WB_0.5.db))
- A file with the covariance matrices of the SNPs within each gene model (such as [this one](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/covariance.DGN-WB_0.5.txt.gz))
- GWAS results (such as [these](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/GWAS.tar.gz), which are just randomly generated)

You can use precalculated databases, or generate new ones with tools in this repository.
GTEx-based tissues and 1000 Genomes covariances precalculated data can be found [here](http://predictdb.hakyimlab.org).
<!-- old box https://app.box.com/s/gujt4m6njqjfqqc9tu0oqgtjvtz9860w  -->
(Please refer to **/software/Readme.md** for more detailed information)

## Setup and Usage Example on a UNIX-like operating system

The following example assumes that you have **python 2.7**, **numpy**, and **scipy** installed.

1) Clone this repository.
```bash
$ git clone https://github.com/hakyimlab/MetaXcan
~~```~~

2) Go to the software folder.
```bash
$ cd MetaXcan/software
```

3) Download sample [data](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/sample_data.tar.gz):
```bash
# You can click on the link above or type the following at a terminal
$ wget https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/sample_data.tar.gz
```
This may take a few minutes depending on your connection: it has to download approximately 200Mb worth of data.
Downloaded data will include an appropiate **Transcriptome Model Database**, a **GWAS/Meta Analysis summary statistics**, and **SNP covariance matrices**.

Extract it with:
```bash
tar -xzvpf sample_data.tar.gz
```

4) Run the High-Level MetaXcan Script
```bash
$ ./MetaXcan.py \
--beta_folder intermediate/beta \
--weight_db_path data/DGN-WB_0.5.db \
--covariance data/covariance.DGN-WB_0.5.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--compressed \
--beta_column BETA \
--pvalue_column P \
--output_file results/test.csv
```
This should take less than a minute on a 3GHZ computer.
Bear in mind that this will generate intermediate data at `intermediate/beta`.
This folder's content's are reused on different runs, not deleted:
you might want to delete this folder before running MetaXcan again,
or specify a different folder on each run.

The example command parameters mean:

* *--beta_folder* Folder where intermediate statistics from the GWAS files will be written to.
* *--weight_db_path* Path to tissue transriptome model
* *--covariance* Path to file containing covariance information. This covariance should have information related to the tissue transcriptome model.
* *--gwas_folder* Folder containing GWAS summary statistics data.
* *--gwas_file_pattern* This option allows the program to select which files from the input to use based on their name.
...This allows to ignore several support files that might be generated at your GWAS analysis, such as plink logs.
* *--beta_column* Tells the program the name of a column containing -phenotype beta data for each SNP- in the input GWAS files.
* *--pvalue_column* Tells the program the name of a column containing -PValue for each SNP- in the input GWAS files.
* *--compressed* This options tells that the input files are in gzip compressed form.
* *--output_file* Path where results will be saved to.

Its output is a CSV file that looks like:

```
gene,gene_name,zscore,effect_size,pvalue,VAR_g,pred_perf_R2,pred_perf_p,pred_perf_q,n_snps_used,n_snps_in_cov,n_snps_in_model
ENSG00000150938,CRIM1,-4.19069760413,-0.231471478373,2.78098095601e-05,0.0983344808163,0.13320775358,1.97496173512e-30,7.47907447189e-30,37,37,37
...
```
Where each row is a gene's association result:
* gene: a gene's id: as listed in the Tissue Transcriptome model.
Ensemble Id for some, while some others (mainly DGN Whole Blood) provide [Genquant](http://www.gencodegenes.org/)'s gene name
* gene_name: gene name as listed by the Transcriptome Model, generally extracted from Genquant
* zscore: MetaXcan's association result for the gene
* effect_size: MetaXcan's association effect size for the gene
* pvalue: P-value of the aforementioned statistic.
* pred_perf_r2: R2 of tissue model's correlation to gene's measured transcriptome (prediction performance)
* pred_perf_pval: pval of tissue model's correlation to gene's measured transcriptome (prediction performance)
* pred_perf_qval: qval of tissue model's correlation to gene's measured transcriptome (prediction performance)
* n_snps_used: number of snps from GWAS that got used in MetaXcan analysis
* n_snps_in_cov: number of snps in the covariance matrix
* n_snps_in_model: number of snps in the model
* var_g: variance of the gene expression, calculated as *W' * G * W*
(where *W* is the vector of SNP weights in a gene's model,
*W'* is its transpose, and *G* is the covariance matrix)

MetaXcan supports a large amount of command line parameters.
Check the Github's ' wiki for those that work best for your data,
and interpreting the results.

## MetaXcan on windows

Please see the following [article](https://github.com/hakyimlab/MetaXcan/wiki/MetaXcan-on-Windows) in the wiki.

## Installation

You also have the option of installing the MetaXcan package to your python distribution.
This will make the **metax** library available for development, and install on your system path
the main MetaXcan scripts.

You can install it from the **software** folder with:

```bash
# ordinary install
$ python setup.py install
```

Alternatively, if you are going to modify the sources, the following may be more convenient:

```bash
# developer mode instalation
python setup.py develop
```

PIP support coming soon-ish.

## Support & Community

Issues and questions can be raised at this repository's issue tracker.

There is also a [Google Group](https://groups.google.com/forum/?hl=en#!forum/predixcanmetaxcan) mail list for general discussion, feature requests, etc. 
Join if you want to be notified of new releases, feature sets and important news concerning this software.

### Cautionary Warning to Existing Users on Updates and Transcriptome Models

Transcriptome Models are a key component of MetaXcan (and PrediXcan) input. As models are improved,
sometimes the format of this databases need be changed too. We only provide support for the very latest databases;
if a user updates their repository to the latest version and MetaXcan complains about the transcriptome weight dbs,
please check if new databases [have been published here](http://predictdb.hakyimlab.org).

For the time being, the only way to use old transcriptome models is to use older versions of MetaXcan.

## Where to go from here

Check [this](https://github.com/hakyimlab/MetaXcan/tree/master/software) if you want to learn more
about more general or advanced usages of MetaXcan.

Check out the [Wiki](https://github.com/hakyimlab/MetaXcan/wiki) for exhaustive usage information.

You will find the manuscript with the theory and rationale for the method at
```bash
/manuscript
```

The code lies at
```bash
/software
```



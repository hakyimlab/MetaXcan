
# MetaXcan

MetaXcan is a set of tools to perform integrative gene mapping studies.
Almost all of the software here is command-line based.

This software has been recently migrated to **python 3** as **python 2** has been sunset.

## S-PrediXcan

S-PrediXcan is an extension of [PrediXcan](https://github.com/hakyimlab/PrediXcan), that infers PrediXcan's results using only summary statistics. It is a component of MetaXcan.
A manuscript desribing S-PrediXcan and the MetaXcan framework and application can be found [here](http://www.biorxiv.org/content/early/2017/05/21/045260).

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

To run S-PrediXcan, you need  [Python 3.5](https://www.python.org/) or higher, with the following libraries:
* [numpy (>=1.14.2)](http://www.numpy.org/)
* [scipy (>=1.2.2)](http://www.scipy.org/) 
* [pandas (>=0.22.0)](http://pandas.pydata.org/)
* [mock](https://github.com/testing-cabal/mock) and [sqlalchemy](https://www.sqlalchemy.org/) are needed for the unit tests.

To run PrediXcan and MulTiPrediXcan, you also need:
* [patsy (>=0.5.0)](https://patsy.readthedocs.io/en/latest/)
* [statsmodels (>=0.10.0)](https://www.statsmodels.org/stable/index.html)
* [h5py (>=2.7.1)](https://github.com/h5py/h5py)
* [h5py-cache (>=1.0.0)](https://pypi.python.org/pypi/h5py-cache/1.0)

[R](https://www.r-project.org/) with [ggplot](http://ggplot2.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) 
is needed for some optional statistics and charts.

# Overview

S-PrediXcan is concerned with obtaining gene-level association tests from ordinary GWAS data.

Ordinarily, a user would need to obtain/download support data sets comprising of:
- A Transcriptome Prediction Model database (an example is [in the example data](https://uchicago.box.com/s/us7qhue3juubq66tktpogeansahxszg9))
- A file with the covariance matrices of the SNPs within each gene model (such as in [the example data](https://uchicago.box.com/s/us7qhue3juubq66tktpogeansahxszg9))

And use them to run S-PrediXcan analysis on:
- GWAS results (such as in  [the example the data](https://uchicago.box.com/s/us7qhue3juubq66tktpogeansahxszg9), from  a randomly generated phenotype)

However, if you have access to interesting data,
you can build your own Transcriptome Prediction Model database 
(using tools such as [PredictDB](http://predictdb.org)), 
and your own covariance (using tools from this repository).

More detailed information and usage documentation can be found at the [Wiki](https://github.com/hakyimlab/MetaXcan/wiki),
where more specific coverage is given to the tools contained here.

# Brief Tour

Scripts **M00_prerequisites.py**, **M01_covariances_correlations.py**, **M02_variances.py**, **M03_betas.py**
and **M04_zscores.py** are steps in a MetaXcan pipeline, although tipically a user will only need to
use the last two. For ease of use, there is a script that performs these last two steps at once, **MetaXcan.py**.

The following scripts are considereded deprecated, 
and exist only for reference purposes:
- **MetaXcanUI.py**, that offers a simplified operation for **MetaXcan.py**'s functionality.
- **M00_prerequisites.py**: Filter down individual dosage data, and convert to supported 'PrediXcan' format.
- **M01_covariances_correlations.py**: builds LD references out of the products from `MOO_prerequisites.py`

All of these runable scripts take different number of command line parameters. 
Run them with **--help** option to see the options, or check the [Wiki](https://github.com/hakyimlab/MetaXcan/wiki) for examples.


## M03_betas.py

This is an optional script that takes GWAS results and converts them into a format that **M04_zscores.py**
can understand, and filters down snps by removing those excluded from the genetic expression.
If you use tools such as [Plink](https://www.cog-genomics.org/plink2), you will need to do this.

This supports both compressed and uncompressed text files as input. When running this script,
you need to specify the column names of expected data. We recommend files that include SNP's GWAS results for
(*beta* or *odd ratio* or *sign of beta*) and p-value. So that, for example, you provide **-or_column** and **--pvalue_column**
when running it.

## M04_zscores.py

This performs the actual S-PrediXcan calculation. It needs a genetic expression model database,
the covariance matrices from this model measured on a specific population, and GWAS data
as in the output of **M03_betas.py**.

If you have the model database and covariance matrix, you can just call this at the end
of your gwas pipeline (optionally using **M03_betas.py**).

## MetaXcan.py

Thin wrapper over **M03_betas.py** and **M04_zscores.py**; it is merely a convenience for running the most frequent S-PrediXcan's steps in a single command.

## MetaXcanUI.py

This script launches a desktop app that processes GWASinput and produces MetaXcan association results.
It is basically a friendlier way to use functionality at **M03_betas.py** and **M04_zscores.py**.
It is deprecated and doesn't support all of the options in the previous scripts (that as a matter of fact haven't been covered in this readme)
but might be friendlier to non technical users.

## MetaMany.py

MetaMany is a script that serially performs S-PrediXcan analysis on a GWAS data set using on multiple tissue models in a single command.

## PrediXcan.py

This is a minimalistic implementation of PrediXcan method. It takes gene expression values (from a plain-text file or HDF5) and computes association to a phenotype vector. 
It supports naive correction for arbitrary covariates (taken from a plain text file or HDF5). there are tools to generate these files [in the main PrediXcan repository](https://github.com/hakyim/PrediXcan).

## MulTiXcan.py

This is the main implementation of the Multi-Tissue association method. 
It takes a collection of files containing gene expression (each file assumed to be a single tissue/study) and computes joint association to a vector phenotype. [The main PrediXcan repository](https://github.com/hakyim/PrediXcan) contaisn tools to generate gene expression files as used by this script.

## SMulTiXcan.py

This is the summary-data based implementation of MulTiXcan. 
You need to run `MetaXcan.py` first for the collection of tissues you are interested in, and a reference LD file such as [this one for v6p models](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/multixcan/snp_covariance.txt.gz)

## Useful Data & Prediction models

We make available several transcriptome predictione models and LD references [here](http://predictdb.org).
These files should be enough for running **MetaXcan.py**, **MulTiXcan.py** and **SMulTiXcan.py** on practically any GWAS study.
we provide a end-to-end [tutorial](https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS), 
for integrating GWAS summary statistics on the latest release of GTEx models.

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


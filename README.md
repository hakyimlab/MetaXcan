> NEW VERSION: we are now correcting the inflation in TWAS due to the pervasive polygenicity of most complex traits. You will need to download the new release of the summary predixcan software as well as the new prediction models, which include necessary correction factors.
> See [TWAS inflation paper here](https://doi.org/10.1101/2023.10.17.562831)

# MetaXcan

MetaXcan is a set of tools to integrate genomic information of biological mechanisms with complex traits.
Almost all of the software here is command-line based.

This software has been recently migrated to **python 3** as **python 2** has been sunset.

## Tools

Here you can find the latest implementation of [PrediXcan]([PrediXcan](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552594/)): **PrediXcan.py**.
This uses individual-level genotype and phenotype, along a mechanism's prediction model (e.g. models predicting expression or splicing quantification),
to compute associations between omic features and a complex trait.

S-PrediXcan is an extension that infers PrediXcan's results using only summary statistics, implemented in **SPrediXcan.py**.
A manuscript describing S-PrediXcan and the MetaXcan framework with an application can be found  [S-PrediXcan](https://www.ncbi.nlm.nih.gov/pubmed/29739930).

[MultiXcan](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007889) (**MulTiXcan.py**) and S-MultiXcan
(**SMulTiXcan.py**) compute omic associations, integrating measurements across tissues while factoring correlation.
For example, if you have prediction models, each trained on different regions of the brain,
MulTiXcan will combine the information across all experiments.
This is effectively a meta-analysis across tissues, where each tissue is an experiment and we explictly account for correlation.


## Prerequisites

The software is developed and tested in Linux and Mac OS environments. The main S-PrediXcan script is also supported in Windows.

To run S-PrediXcan, you need  [Python 3.5](https://www.python.org/) or higher, with the following libraries:
* [numpy (>=1.11.1)](http://www.numpy.org/)
* [scipy (>=0.18.1)](http://www.scipy.org/) 
* [pandas (>=0.18.1)](http://pandas.pydata.org/)
* [sqlalchemy](https://www.sqlalchemy.org/) is needed at some unit tests.

To run PrediXcan Associations and MulTiXcan, you also need:
* [patsy (>=0.5.0)](https://patsy.readthedocs.io/en/latest/)
* [statsmodels (>=0.8.0)](https://www.statsmodels.org/stable/index.html)
* [h5py (>=2.7.1)](https://github.com/h5py/h5py)
* [h5py-cache (>=1.0.0)](https://pypi.python.org/pypi/h5py-cache/1.0) *Now folded into h5py

To run prediction of biological mechanisms on individual-level data, you will also need:
* [bgen_reader (>=3.0.3)](https://pypi.org/project/bgen-reader/)
* [cyvcf2 (>=0.8.0)](https://brentp.github.io/cyvcf2)

[R](https://www.r-project.org/) with [ggplot](http://ggplot2.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) 
is needed for some optional statistics and charts.

We recommend a tool like [Conda](https://www.anaconda.com/distribution/) to set up a working environment for MetaXcan.
Tools like [pyenv](https://github.com/pyenv/pyenv) also work, but the **bgen-reader** dependency currently takes some effort to get going on **pyenv**.

## Example conda environment setup

A quick-and-dirty solution to install the basic requirements is using [Miniconda](https://www.anaconda.com/open-source) and the file [`software/conda_env.yaml`](https://github.com/hakyimlab/MetaXcan/blob/master/software/conda_env.yaml)
in this repository to create a working environment.

```bash
conda env create -f /path/to/this/repo/software/conda_env.yaml
conda activate imlabtools
```

## Useful Data & Prediction models

We make available several transcriptome predictione models and LD references [here](http://predictdb.org).
These files should be enough for running **SPrediXcan.py**, **MulTiXcan.py** and **SMulTiXcan.py** on practically any GWAS study.
We highly recommend **MASHR** models therein, as they are parsimonious and biologically-informed, using fine-mapped variants and cross-tissue QTL patterns.
In the following we use `gene` recurrently to refer to the prediction model of a genetic feature,
but it can stand for other units such as prediction of an intron's quantification.

we provide a end-to-end [tutorial](https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS), 
for integrating GWAS summary statistics on the latest release of GTEx models.

## Project Layout

**software** folder contains an implementation of S-PrediXcan's method and associated tools. 
The following scripts from that folder constitute different components in the MetaXcan pipeline:

```bash
SPrediXcan.py
PrediXcan.py
MulTiXcan.py
SMulTiXcan.py
```
, although `SPrediXcan.py` is the most widely applicable. `SPrediXcan.py` script contains the current implementation of S-PrediXcan. 
`MulTiXcan.py` and `SMulTiXcan.py` are the multiple-tissue methods.
`MultiXcan.py` uses as input the predicted levels generated by `PrediXcan.py`.

The rest of the scripts in **software** folder are python packaging support scripts,
and convenience wrappers such as the GUI.

Subfolder **software/metax** contains the bulk of Metaxcan's logic, implemented as a python module.

## S-PrediXcan Input data

S-PrediXcan will calculate the gene-level association results from GWAS summary statistics. 
It supports most GWAS formats by accepting command line argument specifying data columns.
Some precalculated data is needed, that must be set up prior to S-PrediXcan execution.

The gist of S-PrediXcan's input is:
- A Transcriptome Prediction Model database (an example is [here](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/DGN-WB_0.5.db))
- A file with the covariance matrices of the SNPs within each gene model (such as [this one](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/covariance.DGN-WB_0.5.txt.gz))
- GWAS results (such as [these](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/GWAS.tar.gz), which were computed on a randomly generated phenotype). 
GWAS results can belong to a single file or be split into multiple ones (i.e. split by chromosome).
You can specify the necessary columns via command line arguments (i.e. which column holds snps, which holds p-values, etc)

You can use precalculated databases, or generate new ones with tools available in [PredictDB repository](https://github.com/hakyimlab/PredictDBPipeline).
GTEx-based tissues and 1000 Genomes covariances precalculated data can be found [here](http://predictdb.org).
<!-- old box https://app.box.com/s/gujt4m6njqjfqqc9tu0oqgtjvtz9860w  -->
(Please refer to **/software/Readme.md** for more detailed information)

###  GWAS summary statistic format

S-PrediXcan  supports a large number of input GWAS formats through command line arguments. By specifying the appropriate 
input file column name, S-PrediXcan will analize the file without extra need for input conversion. Input GWAS files can be plain text files or gzip-compressed.

For example, you can specify an effect allele column and a standard error column, or a pvalue column and an odds ratio column, or only a GWAS zscore column. 
S-PrediXcan will try  to use the following (in that order) if available from the command line arguments and input GWAS file:
 
1) use a z-score column if available from the arguments and input file;
2) use a p-value column and either effect, odd ratio or direction column;
3) use effect size (or odd ratio) and standard error columns if available.

Check the Github's ' wiki for those that work best for your data,
and interpreting the results. For example, if your GWAS has p-values that are too small (i.e 1e-350),
then you should avoid specifying a p-value column because numerical problems might arise; you should use effect size and standard error instead.

## A remark on individual-level genotype format

PrediXcan supports three input file formats:
- vcf
- bgen
- internal "dosage format".

Associations are output as a tab-separated file.

Predicted levels can be output as both text files or HDF5 files.
HDF5 files allow a more efficient computation of MultiXcan,
as only data for a single gene/inton/whaever across all tissues can be loaded at a time.


## Setup and Usage Example on a UNIX-like operating system

The following example assumes that you have **python 3.5** (or higher), **numpy**, and **scipy** installed.

1) Clone this repository.
```bash
$ git clone https://github.com/hakyimlab/MetaXcan
```

2) Go to the software folder.
```bash
$ cd MetaXcan/software
```

3) Download example [data](https://uchicago.box.com/s/us7qhue3juubq66tktpogeansahxszg9).

This may take a few minutes depending on your connection: it has to download approximately 200Mb worth of data.
Downloaded data will include an appropiate **Transcriptome Model Database**, a **GWAS/Meta Analysis summary statistics**, and **SNP covariance matrices**.

Extract it with:
```bash
tar -xzvpf sample_data.tar.gz
```

4) Run the High-Level S-PrediXcan Script
```bash
./SPrediXcan.py \
--model_db_path data/DGN-WB_0.5.db \
--covariance data/covariance.DGN-WB_0.5.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/test.csv
```
This should take less than a minute on a 3GHZ computer. For the full specification of command line parameters, you can check the [wiki](https://github.com/hakyimlab/MetaXcan/wiki/MetaXcan's-Command-Line-Reference).


The example command parameters mean:

* `--model_db_path` Path to tissue transriptome model
* `--covariance` Path to file containing covariance information. This covariance should have information related to the tissue transcriptome model.
* `--gwas_folder` Folder containing GWAS summary statistics data.
* `--gwas_file_pattern` This option allows the program to select which files from the input to use based on their name.
...This allows to ignore several support files that might be generated at your GWAS analysis, such as plink logs.
* `--snp_column` Argument with the name of the column containing the RSIDs.
* `--effect_allele_column` Argument with the name of the column containing the effect allele (i.e. the one being regressed on).
* `--non_effect_allele_column` Argument with the name of the column containing the non effect allele.
* `--beta_column` Tells the program the name of a column containing -phenotype beta data for each SNP- in the input GWAS files.
* `--pvalue_column` Tells the program the name of a column containing -PValue for each SNP- in the input GWAS files.
* `--output_file` Path where results will be saved to.

Its output is a CSV file that looks like:

```
gene,gene_name,zscore,effect_size,pvalue,var_g,pred_perf_r2,pred_perf_pval,pred_perf_qval,n_snps_used,n_snps_in_cov,n_snps_in_model
ENSG00000150938,CRIM1,4.190697619877402,0.7381499095142079,2.7809807629839122e-05,0.09833448081630237,0.13320775358,1.97496173512e-30,7.47907447189e-30,37,37,37
...
```
Where each row is a gene's association result:
* `gene`: a gene's id: as listed in the Tissue Transcriptome model.
Ensemble Id for most gene model releases. Can also be a intron's id for splicing model releases.
* `gene_name`: gene name as listed by the Transcriptome Model, typically HUGO for a gene. It can also be an intron's id.
* `zscore`: S-PrediXcan's association result for the gene, typically HUGO for a gene.
* `effect_size`: S-PrediXcan's association effect size for the gene. Can only be computed when `beta` from the GWAS is used.
* `pvalue`: P-value of the aforementioned statistic.
* `pred_perf_r2`: (cross-validated) R2 of tissue model's correlation to gene's measured transcriptome (prediction performance). Not all model families have this (e.g. MASHR).
* `pred_perf_pval`: pval of tissue model's correlation to gene's measured transcriptome (prediction performance). Not all model families have this (e.g. MASHR).
* `pred_perf_qval`: qval of tissue model's correlation to gene's measured transcriptome (prediction performance). Not all model families have this (e.g. MASHR).
* `n_snps_used`: number of snps from GWAS that got used in S-PrediXcan analysis
* `n_snps_in_cov`: number of snps in the covariance matrix
* `n_snps_in_model`: number of snps in the model
* `var_g`: variance of the gene expression, calculated as `W' * G * W`
(where `W` is the vector of SNP weights in a gene's model,
`W'` is its transpose, and `G` is the covariance matrix)

If `--additional_output` is used when running S-PrediXcan, you'll get two additional columns:
* `best_gwas_p`: the highest p-value from GWAS snps used in this model
* `largest_weight`: the largest (absolute value) weight in this model

## S-PrediXcan on windows

Please see the following [article](https://github.com/hakyimlab/MetaXcan/wiki/MetaXcan-on-Windows) in the wiki.

## S-MultiXcan

Run `SMultiXcan.py --help` to see arguments and options.
Output is a tab-separated text file with the following columns:
* `gene`: a gene's id: as listed in the Tissue Transcriptome model.
Ensemble Id for most gene model releases. Can also be a intron's id for splicing model releases.
* `gene_name`: gene name as listed by the Transcriptome Model, typically HUGO for a gene. It can also be an intron's id.
* `pvalue`: significance p-value of S-MultiXcan association
* `n`: number of "tissues" available for this gene
* `n_indep`: number of independent components of variation kept among the tissues' predictions. (Synthetic independent tissues)
* `p_i_best`: best p-value of single-tissue S-PrediXcan association.
* `t_i_best`: name of best single-tissue S-PrediXcan association.
* `p_i_worst`: worst p-value of single-tissue S-PrediXcan association.
* `t_i_worst`: name of worst single-tissue S-PrediXcan association.
* `eigen_max`: In the SVD decomposition of predicted expression correlation: eigenvalue (variance explained) of the top independent component
* `eigen_min`: In the SVD decomposition of predicted expression correlation: eigenvalue (variance explained) of the last independent component
* `eigen_min_kept`: In the SVD decomposition of predicted expression correlation: eigenvalue (variance explained) of the smalles independent component that was kept.
* `z_min`: minimum z-score among single-tissue S-Predican associations.
* `z_max`: maximum z-score among single-tissue S-Predican associations.
* `z_mean`: mean z-score among single-tissue S-Predican associations.
* `z_sd`: standard deviation of the mean z-score among single-tissue S-Predican associations.
* `tmi`: trace of `T * T'`, 
where `T`is correlation of predicted expression levels for different tissues 
multiplied by its SVD pseudo-inverse. 
It is an estimate for number of indepent components of variation in predicted expresison across tissues (typically close to `n_indep`)
* `status`: If there was any error in the computation, it is stated here


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

You can check [here](https://github.com/hakyimlab/MetaXcan/wiki) for the release history.

### Cautionary Warning to Existing Users on Updates and Transcriptome Models

Transcriptome Models are a key component of PrediXcan and S-PrediXcan input. As models are improved,
sometimes the format of these databases needs be changed too. We only provide support for the very latest databases;
if a user updates their repository clone to the latest version and MetaXcan complains about the transcriptome weight dbs,
please check if new databases [have been published here](http://predictdb.org).

For the time being, the only way to use old transcriptome models is to use older versions of MetaXcan.

Older versions of MetaXcan have a **MetaXcan.py** script, when it meant to be an entry point to all MetaXcan tools, but it has since been renamed **SPrediXcan.py**.

## Where to go from here

Check [software](https://github.com/hakyimlab/MetaXcan/tree/master/software) folder in this repository if you want to learn
about more general or advanced usages of S-PrediXcan, or MulTiXcan and SMulTiXcan.

Check out the [Wiki](https://github.com/hakyimlab/MetaXcan/wiki) for exhaustive usage information.

The code lies at
```bash
/software
```

New release and features coming soon!



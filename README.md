
# msqrob2TMT paper

This repository contains the code used to generate the results and
generate the manuscript published in:

> Vandenbulcke, S., Vanderaa, C., Crook, O., Martens, L. & Clement, L.
> msqrob2TMT: robust linear mixed models for inferring differential
> abundant proteins in labelled experiments with arbitrarily complex
> design. bioRxiv 2024.03.29.587218 (2024)
> doi:10.1101/2024.03.29.587218

We also provide 2 tutorial vignettes that illustrate how to apply
msqrob2TMT in practice.

## Abstract

Labelling strategies in mass spectrometry (MS)-based proteomics
enhance sample throughput by enabling the acquisition of multiplexed
samples within a single run. However, contemporary experiments often
involve increasingly complex designs, where the number of samples
exceeds the capacity of a single run, resulting in a complex
correlation structure that must be addressed for accurate statistical
inference and reliable biomarker discovery. To this end, we introduce
msqrob2TMT, a suite of mixed model-based workflows specifically
designed for differential abundance analysis in labelled MS-based
proteomics data. msqrob2TMT accommodates both sample-specific and
feature-specific (e.g., peptide or protein) covariates, facilitating
inference in experiments with arbitrarily complex designs and allowing
for explicit correction of feature-specific covariates. We benchmark
our innovative workflows against state-of-the-art tools, including
DEqMS, MSstatsTMT, and msTrawler, using two spike-in studies. Our
findings demonstrate that msqrob2TMT offers greater flexibility,
improved modularity, and enhanced performance, particularly through
the application of robust ridge regression. Finally, we demonstrate
the practical relevance of msqrob2TMT in a real mouse study,
highlighting its capacity to effectively account for the complex
correlation structure in the data.

## Datasets

To demonstrate the performance and application of msqrob2TMT, we used
3 datasets.

### spikein1

The spikein1 dataset has been published by [Huang et al.
2020](http://dx.doi.org/10.1074/mcp.RA120.002105). It has the
following design: UPS1 peptides at concentrations of 500, 333, 250,
and 62.5 fmol were spiked into 50 µg of SILAC HeLa peptides. This
series forms a dilution gradient of 1, 0.667, 0.5, and 0.125 relative
to the highest UPS1 peptide concentration (500 fmol). A reference
sample was prepared by combining the diluted UPS1 peptide samples with
50 µg of SILAC HeLa peptides. Each dilution and the reference sample
were processed in duplicate, yielding a total of ten samples. These
samples were subsequently labelled using TMT10-plex reagents and
combined for LC-MS/MS analysis, hereafter referred to as a "mixture."
This protocol was repeated five times. Each mixture was analysed in
technical triplicate, resulting in a total of 15 MS runs. This
experimental design simulates a scenario with 10 biological replicates
per condition, consisting of 2 replicates per condition within each of
the five mixtures.

### spikein2

The spikein1 dataset has been published by [O'Brien et al.
2024](http://dx.doi.org/10.1038/s41592-023-02120-6). It is an
experiment which consists of a constant mouse plasma background in
which yeast proteins were then spike-in with the following dilution
ratios 1, 2, 3, 5, 11, 14, 18, 20, 24, 28 and 32. The yeast proteins
were harvested from different media; glycerol + ammonium sulfate,
galactose + urea, galactose + monosodium glutamate, galactose +
ammonium sulfate, glucose + ammonium sulfate, glucose + monosodium
glutamate, and glucose + urea. The samples were labelled using the
TMTpro 18-plex and pooled in six mixtures, each containing 15 samples.
Two of the remaining TMT channels are reference channels, one
reference channel is the combination of one batch and the other
reference channel consists of all batches. For this dataset the ground
truth is known: a true positive is a yeast protein that is returned as
significant while a significant mouse protein is a false positive.

### mouse

The mouse dataset has been published by [Plubell et al.
2017](https://doi.org/10.1074/mcp.m116.065524). Twenty mice were
divided into groups to investigate the effects of low-fat (LF) and
high-fat (HF) diets. The experiment involves two factors: the type of
diet (LF or HF), and the duration of the diet, categorised as either
short-term (8 weeks) or long-term (18 weeks). This design resulted in
four distinct groups corresponding to each combination of diet type
and duration. Five mice, i.e. biological replicates, were randomly
assigned to each group and samples from their epididymal adipose
tissue were randomly assigned to three TMT10-plex mixtures. Within
each mixture, two reference channels were included, each containing
pooled samples representing a range of peptides from all samples. Not
all TMT channels were utilised, leading to an unbalanced experimental
design. Finally, each TMT mixture was fractionated into nine parts and
analysed using synchronous precursor selection, culminating in a total
of 27 MS runs.

## Repository organisation

The repository contains 4 main folders:

### `vignettes/`

To illustrate how to use msqrob2TMT, we provide 2 tutorial vignettes:

- `spikein1.qmd`: the vignette presents the theoretical aspects and the
  practical implementation of labelling-based MS proteomics data
  analysis side-by-side. As this tutorial focuses on a spike-in
  dataset, the impact of modelling decision on the accuracy of the
  method are illustrated.
- `mouse.qmd`: the vignette illustrates the application of msqrob2TMT
  on a real-life dataset, showing how to translate biological
  questions into statistical model
  definition.

### `scripts/`

This folder contains all the code used to preprocess the data, run the
different modelling methods and to generate the figures from the
manuscript. Scripts have been organised by dataset, one folder for
each dataset.

### `figs/`

This folder contains all the figures generated by the scripts and
included in the manuscript.

### `data/`

This folder will contain intermediate tables generated by the scripts.
This folder is empty on GitHub, but can be populated with the files in
the `processed.zip` archive provided on our [Zenodo
repository](https://zenodo.org/records/14767905).

## Getting the data

Although all datasets are publicly available, we have deposited all
the required data to run the analyses in a
[Zenodo repository](https://zenodo.org/records/14767905).

All input data are automatically downloaded within the scripts, using
the Bioconductor `BiocFileCache` package. This package will ensure an
efficient file management, downloading the data only once without the
need for user intervention.

## Setting up the computational environment

There are two strategies for obtaining the computational environment
(i.e. the set of installed software) to reproduce our findings: using
Docker (recommended) or going for a manual setup. For both strategies,
you must first clone this GitHub repository on your machine either by
downloading the code as a zip file or through `git`.

#### Docker

You can reproduce the paper's figures on your local machine using
`Docker`. You must first install Docker, then pull the image from the
[DockerHub
repository](https://registry.hub.docker.com/r/cvanderaa/msqrob2tmt).

```
docker pull cvanderaa/msqrob2tmt
```

Then, you can start an Rstudio session within a Docker container using:

```
docker run \
    -e PASSWORD=bioc \
    -p 8787:8787 \
    -v `pwd`:/home/rstudio/ \
    cvanderaa/msqrob2tmt
```

Note you should use `%CD%` instead of `pwd` when using Windows.

Open your browser and go to http://localhost:8787. The USER is
`rstudio` and the password is `bioc`. See the [DockerHub
repository](https://hub.docker.com/repository/docker/cvanderaa/msqrob2tmt)
for more detailed information on getting started with `Docker`.

Similarly, you can run the Docker environment from the command
line, providing more flexibility for developers and data
analysts.

```
docker run -v `pwd`:/home/rstudio/ -it cvanderaa/msqrob2tmt bash
```

You can replace `bash` by your favourite shell program or by `R`
to immediately start an R environment.

#### Note to future self: build the docker image

If new dependencies are required, update the `Dockerfile` accordingly.
Then build the image (make sure to `cd` in the repo):

```
docker build -t cvanderaa/msqrob2tmt .
```

See the [DockerHub
repo](https://registry.hub.docker.com/r/cvanderaa/msqrob2tmt) for more
information. When complete, push the new image to DockerHub:

```
docker push cvanderaa/msqrob2tmt
```

#### Manual setup

Alternatively, you can set the software environment yourself. To
reproduce the article and its results, you'll need the following
prerequisites:

- Install `R` (suggest R version >= 4.4.0)
- Install Bioconductor using

```r
install.packages("BiocManager")
```

Make sure that `BiocManager::version()` returns at least 3.20.

- Install [Quarto](https://quarto.org/)
- Install the following R packages, make sure `msqrob2` is version >=
  1.10:

```r
BiocManager::install(c("BiocParallel", "BiocStyle", "Biostrings", "curl", "DEqMS", "dplyr", "GOfuncR", "QFeatures", "lme4", "msqrob2", "MSstatsTMT", "Mus.musculus", "pryr", "quarto", "readr", "stringr", "tidyverse", "cowplot", "ggVennDiagram", "ggridges", "gt", "flextable", "ggh4x", "gprofiler2", "UpSetR", "ggborderline", "viridisLite", "impute", "ExploreModelMatrix"))
```

Make sure that `msqrob2` is at least version 1.14.1 and `QFeatures` is
at least version 1.16.1. If this is not the case, you can install the
latest version using

```r
BiocManager::install(c(
  "rformassspectrometry/QFeatures",
  "statsOmics/msqrob2"
))
```

Install `msTrawler` time stamped at commit `078eb79`.

```r
BiocManager::install("calico/msTrawler", ref = "078eb79")
```

You should be ready to reproduce the results, although this approach
can be error-prone unlike the Docker approaches.

## Running the tutorials

Once you have setup the computational environment, you can run our
tutorial vignettes for the spikin1 and the mouse dataset.

## Reproducing the results

Once you have setup the computational environment, you can also rerun
our analyses. Each dataset can be analysed independently using the
scripts from the corresponding subfolder in`scripts/`. The scripts
will generate intermediate data files stored in `data/` that serve as
input for following scripts. Therefore, scripts should be run in the
order shown below. However, the intermediate files stored in `data/`
are also provided on our [Zenodo
repository](https://zenodo.org/records/14767905). If you download the
`processed.zip` and unzip it in `data/`, the order no longer matters.

#### spikein1

- spikein1_preprocess.R
- spikein1_deqms.R
- spikein1_msqrob2tmt.R
- spikein1_msstatstmt.R
- spikein1_mstrawler.R
- spikein1_compare_preprocessing.R
- spikein1_figures.R

#### spikein2

- spikein2_preprocess.R
- spikein2_deqms.R
- spikein2_msqrob2tmt.R
- spikein2_msstatstmt.R
- spikein2_mstrawler.R
- spikein2_figures.R

#### mouse

- mouse_msqrob2tmt_mixture.R
- mouse_msqrob2tmt.R
- mouse_msstatstmt.R
- mouse_figures.R

Note that the `*_figures.R` scripts will generate figures in the
`figs/` folder.

## Citation

> Vandenbulcke, S., Vanderaa, C., Crook, O., Martens, L. & Clement, L.
> msqrob2TMT: robust linear mixed models for inferring differential
> abundant proteins in labelled experiments with arbitrarily complex
> design. bioRxiv 2024.03.29.587218 (2024)
> doi:10.1101/2024.03.29.587218

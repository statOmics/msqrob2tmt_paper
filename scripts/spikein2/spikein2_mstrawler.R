####---- msTrawler workflow analyses on the spike-in2 data ----####

## This script runs msTrawler on the spikein2 data from O'Brien et al.
## 2024 where yeast proteins are spiked in a mouse background.

####---- Setup environment ----####

# Read in the relevant libraries
library("msTrawler")
library("BiocFileCache")
library("dplyr")

dataDir <- "data/"

####---- Load data ----####

## Sample annotations
sampleFile <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein2_annotations.csv?download=1") |>
    read.csv()

## Covariate file
covariateFile <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein2_covariateFile.csv?download=1") |>
    read.csv()

spikein2 <- readRDS(paste0(dataDir, "spikein2_input.rds"))
## msTrawler performs imputation with 0.
spikein2[is.na(spikein2)] <- 0

####---- Run msTrawler ----####

curDir <- getwd()
setwd(dataDir)
set.seed(123)
msTrawler::msTrawl(
    spikein2,
    sampleFile = sampleFile,
    covariateFile = covariateFile,
    scaleSN = 1, ## default
    lod = 0.01, ## default
    imputePenalty = 1, ## default
    minAbove = 3, ## as in author's script for multibatch (spikein2) analysis
    ssnFilter = NULL, ## disable, already performed
    outlierCutoff = NULL, ## disable, already performed
    N_SUM = 3, ## default
    swapProtein = FALSE, ## default
    maxPep = 25, ## default
    colAdjust = NULL,## disable, already performed
    colRatios = NULL, ## disable, already performed
    dropContam = FALSE, ## disable, already performed
    dropReverse = FALSE, ## disable, already performed
    peptideAnalysis = FALSE, ## default
    minRE = 5, ## default
    timeDiff = FALSE ## not applicable
)
setwd(curDir)

####---- Retrieve and save results ----####

## Make a table of comparison (used for automating result extraction)
comparisons <- unique(sampleFile$Dilution)[-9] |>
    as.numeric() |>
    sort(decreasing = TRUE) |>
    combn(2) |>
    t()

## Extract statistical inference results from local files (sigh...)
resultFiles <- list.files(dataDir, "^Factor_Dilution.*csv", full.names = TRUE)
inference <- lapply(1:nrow(comparisons), function(i) {
    group1 <- comparisons[i, 1]
    group2 <- comparisons[i, 2]
    out <- paste0("_", group2, ".csv") |>
        grep(resultFiles, value = TRUE) |>
        read.csv()
    colSel <- c(
        "Protein",
        paste0("Est_", group1),
        paste0("Pval_", group1),
        paste0("Qval_", group1)
    )
    out <- out[, colSel]
    colnames(out) <- c("Protein", "logFC", "pval", "adjPval")
    ## Add more columns for consistency with other methods
    out$t <- NA
    out$df <- NA
    out$se <- NA
    out$Comparison <- paste(group1, "-", group2)
    out$Model <- "msTrawler"
    out
})
inference <- do.call(rbind, inference)
inference$Preprocessing <- "msTrawler"

## Store the modelled spikein1 data
saveRDS(inference, paste0(dataDir, "spikein2_model_msTrawler.rds"))

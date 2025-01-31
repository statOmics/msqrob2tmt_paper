
## DEqMS analysis on the spikein2 data ----####

## This script runs DEqMS on the spikein2 data from O'Brien et al.
## 2024 where yeast proteins are spiked in a mouse background.

####---- Setup environment ----####

# Read in the relevant libraries
library("DEqMS")
library("BiocFileCache")
library("tidyr")
library("dplyr")

dataDir <- "data/"

####---- Load data ----####

## Sample annotations
annotations <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein2_annotations.csv?download=1") |>
    read.csv()
annotations$Run <- sub("_.*", "", annotations$SampleID)
annotations$Run <- sub("-", "", annotations$Run)
annotations$SampleID <- sub("Sn", "Adjusted.Intensity", annotations$SampleID)
rownames(annotations) <- annotations$SampleID

## Read PSM data
spikein2 <- readRDS(paste0(dataDir, "spikein2_input.rds")) |>
    ## Drop the S/N normalisation
    select(!contains(".Sn")) |>
    ## Collapse all channels in a single column
    pivot_longer(
        cols = contains("Intensity"),
        names_to = "Channel"
    ) |>
    ## Expand the collapsed column into samples = plex + channel
    pivot_wider(id_cols = c("Protein.ID","Peptide"),
                names_from = c(Plex, Channel),
                values_from = value)

## Remove normalisation channel
normSamples <- annotations$SampleID[annotations$Bridge == 1]
spikein2 <- spikein2[, !colnames(spikein2) %in% normSamples]

## Match PSM data and annotations because one sample is not annotated
quants <- spikein2[, -(1:2)]
quants <- quants[, colnames(quants) %in% annotations$SampleID]
spikein2 <- cbind(spikein2[, 1:2], quants)
annotations <- annotations[colnames(quants), ]

## Get PSM counts per protein for each run, used for spectral count
## scaling by DEqMS
counts <- lapply(split(colnames(spikein2)[-(1:2)], annotations$Run), function(ii) {
    observed <- rowSums(!is.na(spikein2[, ii])) != 0
    proteinsInRun <- spikein2$Protein.ID[observed]
    data.frame(table(Protein.ID = proteinsInRun))
})
counts <- Reduce(
    function(x, y) merge(x, y, by = "Protein.ID", all = TRUE),
    counts
)
countsn <- counts$Protein.ID
counts <- rowMins(as.matrix(counts[, -1]), na.rm = TRUE)
names(counts) <- countsn

####---- Preprocessing ----####

## Log-transform
spikein2 <- mutate_at(spikein2, vars(contains("Intensity")), log2)

## The functions below is copy-pasted from DEqMS::medianSweeping(),
## removing the call to DEqMS::equalMedianNormalization() since
## normalization has already been done beforehand
medianSweepingNoNorm <- function (dat, group_col = 2) {
    dat.ratio = dat
    dat.ratio[, 3:ncol(dat)] = dat.ratio[, 3:ncol(dat)] - matrixStats::rowMedians(as.matrix(dat.ratio[,
                                                                                                      3:ncol(dat)]), na.rm = TRUE)
    dat.summary = plyr::ddply(dat.ratio, colnames(dat)[group_col],
                              function(x) matrixStats::colMedians(as.matrix(x[, 3:ncol(dat)]),
                                                                  na.rm = TRUE))
    colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
    dat.new = dat.summary[, -1]
    rownames(dat.new) = dat.summary[, 1]
    return(dat.new)
}

# spikein2 <- medianSweepingNoNorm(spikein2, group_col = 1)
spikein2 <- lapply(split(annotations$SampleID, annotations$Run), function(cols) {
    x <- spikein2[, c("Peptide", "Protein.ID", cols)]
    x <- medianSweepingNoNorm(x, group_col = 2)
    x$Protein.ID <- rownames(x)
    x
}) |>
    Reduce(f = function(x, y) full_join(x, y, by = join_by(Protein.ID)))
rownames(spikein2) <- spikein2$Protein.ID
spikein2$Protein.ID <- NULL

####---- Run DEqMS ----####

## Make design matrix
## Fixed effect for condition and fixed effect for run
annotations <- annotations[colnames(spikein2), ]
design <- model.matrix(~ 0 + Dilution + Media + Run, data = annotations)

## Make contrasts
paramNames <- grep("Dilution", colnames(design), value = TRUE)
paramValues <- sub("Dilution", "", paramNames) |>
    as.numeric() |>
    sort(decreasing = TRUE)
combinations <- paste0("Dilution", paramValues) |>
    combn(2, paste, collapse = " - ")
contrast <- makeContrasts(contrasts = combinations, levels = design)

## Run DEqMS modelling
fit1 <- lmFit(spikein2, design)
fit2 <- eBayes(contrasts.fit(fit1, contrasts = contrast))
## Note that some proteins have very low df, we will remove them as
## this will lead to data inconsistencies downstream the script
sel <- fit1$df.residual > 1 & rowMeans(is.na(fit2$coefficients)) == 0
fit1 <- lmFit(spikein2[sel, ], design)
fit2 <- eBayes(contrasts.fit(fit1, contrasts = contrast))
fit2$count <- counts[rownames(fit2$coefficients)]
fit3 <- spectraCounteBayes(fit2)

####---- Format and save output ----####

out <- lapply(colnames(fit3$coefficients), function(ii) {
    test <- outputResult(fit3, coef_col = ii)
    ## Spectral count scaling fails, so use unscaled results
    test <- test[, c("logFC", "t", "sca.adj.pval", "sca.P.Value")]
    colnames(test) <- c("logFC", "t", "adjPval", "pval")
    test$Protein <- rownames(test)
    test$df <- NA
    test$se <- NA
    test$Comparison <- gsub("Dilution", "", ii)
    test
})
out <- do.call(rbind, out)
out$Model <- "DEqMS"
out$Preprocessing <- "DEqMS"

## Store the modelled spikein1 data
saveRDS(out, paste0(dataDir, "spikein2_model_DEqMS.rds"))

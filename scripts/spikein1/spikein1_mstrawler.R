####---- msTrawler workflow analyses on the UPS spike-in data ----####

## This script runs msTrawler on the UPS spike-in data.

####---- Setup environment ----####

# Read in the relevant libraries
library("msTrawler")
library("BiocFileCache")
library("dplyr")

dataDir <- "data/"

####---- Get annotations ----####

## msTrawler requires the annotations are provided in a data.frame,
## with 3 columns: Covariate to model, Bridge (0 or 1), and SampleID.
sampleFile <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein1_annotations.csv?download=1") |>
    read.csv()

## msrawler uses File.ID as plex identifier, but this is not present
## in the current sampleFile, so we convert Run (i.e. Spectrum.File)
## in File.ID
file2Id <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein1_psms.txt?download=1") |>
    read.delim() |>
    select(Spectrum.File, File.ID) |>
    unique()
sampleFile$Plex <- file2Id$File.ID[match(sampleFile$Run, file2Id$Spectrum.File)]

## SampleID is the plex + SN normalised column name.
sampleFile$SampleID <- paste0(sampleFile$Plex, "_X", tolower(sampleFile$Channel), ".Sn")
sampleFile$Bridge <- ifelse(sampleFile$Condition == "Norm", 1, 0)
sampleFile <- sampleFile[, c("SampleID", "Condition", "Bridge")]

## Remove one of the two normalisation (bridge) channels, msTrawler
## allows only for one such channel.
sampleFile <- sampleFile[!grepl("131", sampleFile$SampleID), ]

## Generate the covariateFile table
covariateFile <- data.frame(
    Covariate = "Condition",
    Type = "Factor",
    levels = 5,
    timeDegree = 0,
    TimeCategory = 0,
    Circadian  = 0
)

## Make a table of comparison (used for automating result extraction)
comparisons <- unique(sampleFile$Condition)[-1] |>
    sort(decreasing = TRUE) |>
    combn(2) |>
    t()

####---- msTrawler data importer ----####

## `msTrawler::convertPdFile()` is bugged. More specifically, Proteins
## identified with the Mascot Node don't have a
## Master.Protein.Accesssions (encoded as ""). This
## [line](https://github.com/calico/msTrawler/blob/078eb790f33667457feb033521c5e76111f6b75c/R/FileConverter.R#L99)
## removes the entry completely, leading to a shift in column index,
## leading to the protein ID to be removed. This is a serious issue as
## the proteins identified by Mascot are the UPS proteins, hence most
## DA proteins can no longer be identified as DA while they are,
## leading to a drop in performance.

## To included meaningful results for msTrawler, I reimplemented the
## function.
convertPdFileFixed <- function(file) {
    data <- read.delim(
        paste0(dataDir, "spikein1_input_msTrawler.txt"),
        check.names = FALSE
    )
    ## Format the quantitative data
    quantColIndex <- grep("Abundance", colnames(data))
    quants <- data[, quantColIndex]
    quants[is.na(quants)] <- 0
    colnames(quants) <- sub(".*: (.*)", "X\\1.Adjusted.Intensity", tolower(colnames(quants)))
    ## Normalise the data
    snRatio <- data$`Average Reporter S/N`
    noise <- rowMeans(quants) / snRatio
    quantsNorm <- quants / noise
    colnames(quantsNorm) <- sub(".Adjusted.Intensity", ".Sn", colnames(quantsNorm))
    ## Create the table expected by msTrawler
    cleanedPeptide <- gsub("\\[|\\]", "", data$`Annotated Sequence`)
    ## Return data
    data.frame(
        PA.Gene.Symbol = NA,
        Protein.ID = data$`Protein Accessions`,
        Peptide = paste(cleanedPeptide, data$Modifications),
        Plex = data$`File ID`,
        quants,
        quantsNorm
    )
}

####---- msTrawler workflow ----####

## Use either the original or the fixed function
inference <- lapply(c(FALSE, TRUE), function(useFix) {
    ## Read data
    f <- paste0(dataDir, "spikein1_input_msTrawler.txt")
    converter <- ifelse(useFix, convertPdFileFixed, convertPdFile)
    spikein1MsTrawler <- converter(f)

    ## Remove one of the two normalisation (bridge) channels, msTrawler
    ## allows only for one such channel.
    spikein1MsTrawler <- spikein1MsTrawler[, !grepl("131", colnames(spikein1MsTrawler))]
    ## Make sure the sampleFile contains the same SampleID as the psm
    ## data
    sampleFileMatched <- sampleFile[sub("_.*", "", sampleFile$SampleID) %in% unique(spikein1MsTrawler$Plex), ]
    sampleFileMatched <- sampleFileMatched[sub(".*_", "", sampleFileMatched$SampleID) %in% colnames(spikein1MsTrawler), ]

    ## Run modelling
    curDir <- getwd()
    setwd(dataDir)
    set.seed(999)
    msTrawl(
        DF = spikein1MsTrawler,
        sampleFile = sampleFileMatched,
        covariateFile = covariateFile,
        timeDiff = FALSE
    )
    setwd(curDir)

    ## Extract statistical inference results from local files (sigh...)
    resultFiles <- list.files(dataDir, "^Factor_Condition.*csv", full.names = TRUE)
    lapply(1:nrow(comparisons), function(i) {
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
        out$Model <- ifelse(useFix, "msTrawler_fixed", "msTrawler")
        out$technicalReplicates <- FALSE
        out
    }) |> do.call(what = rbind)
})

## Store the modelled spikein1 data
inference <- do.call(rbind, inference)
saveRDS(inference, paste0(dataDir, "spikein1_model_msTrawler.rds"))

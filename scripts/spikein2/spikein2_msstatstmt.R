## MSstatsTMT analysis on the spikein2 data ----####

## This script runs MSstatsTMT on the spikein2 data from O'Brien et al.
## 2024 where yeast proteins are spiked in a mouse background.

####---- Setup environment ----####

# Read in the relevant libraries
library("MSstatsTMT")
library("BiocFileCache")
library("tidyr")
library("dplyr")

dataDir <- "data/"

####---- Load data ----####

## Load sample annotations
annotations <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein2_annotations.csv?download=1") |>
    read.csv()
annotations$SampleID <- sub("Sn", "Adjusted.Intensity", annotations$SampleID)

## Load PSM data
spikein2 <- readRDS(paste0(dataDir, "spikein2_input.rds")) |>
    ## Drop the S/N normalisation
    select(!contains(".Sn")) |>
    ## Collapse all channels in a single column
    pivot_longer(
        cols = contains("Intensity"),
        names_to = "Channel", values_to = "intensity"
    ) |>
    ## Combine with sample annotations
    mutate(SampleID = paste0(Plex, "_", Channel)) |>
    left_join(annotations, join_by(SampleID)) |>
    ## Remove non-annotated samples
    filter(!is.na(Dilution))

## Need to remove _ already in some peptide names, as the msTrawler
## workflow already summed peptides with multiple PSM. Otherwise, this
## leads to an error from a str_split function within
## proteinSummarization
spikein2$Peptide <- gsub("_", "", spikein2$Peptide)

## Remove normalization samples
spikein2$Dilution <- sub("B", "Norm", spikein2$Dilution)

## Convert to a MSstatsTMT-compatible table
spikein2 <- data.frame(
    ProteinName = spikein2$Protein.ID,
    PeptideSequence = spikein2$Peptide,
    Charge = 1,
    PSM = paste0(spikein2$Peptide, "_", spikein2$z),
    Mixture = spikein2$Plex,
    TechRepMixture = 1,
    Run = spikein2$Plex,
    Channel  = spikein2$Channel,
    Condition = spikein2$Dilution,
    BioReplicate = spikein2$Dilution,
    Intensity = spikein2$intensity
)

####---- MSstats modelling ----####

spikein2 <- proteinSummarization(
    spikein2,
    method = "msstats",
    global_norm = FALSE,
    reference_norm = TRUE,
    remove_norm_channel = TRUE,
    remove_empty_channel = TRUE,
    MBimpute = TRUE,
    maxQuantileforCensored = NULL,
    use_log_file = FALSE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL,
    msstats_log_path = NULL
)

####---- Statistical inference ----####

## Make contrast matrix
paramNames <- unique(spikein2$ProteinLevelData$Condition) |>
    as.character() |>
    as.numeric() |>
    sort(decreasing = TRUE)
paramNames <- paste0("Dilution", paramNames)
contrasts <- combn(paramNames, 2, simplify = FALSE) |>
    sapply(function(x) paste(paste(x, collapse = " - "), " = 0")) |>
    msqrob2::makeContrast(paramNames) |>
    t()
dimnames(contrasts) <- lapply(dimnames(contrasts), function(x) gsub("Dilution", "", x))

## Perform statistical inference
tests <- groupComparisonTMT(
    spikein2,
    contrast.matrix = contrasts,
    moderated = TRUE
)$ComparisonResult

####---- Format and save output ----####

tests <- tests[, -8] ## remove the "issue" column
colnames(tests) <- c("Protein", "Comparison", "logFC", "se", "df", "pval", "adjPval")
tests$t <- NA
tests$Model <- "MSstatsTMT"
tests$Preprocessing <- "MSstatsTMT"

## Store the modelled spikein1 data
saveRDS(data.frame(tests), paste0(dataDir, "spikein2_model_MsstatsTMT.rds"))

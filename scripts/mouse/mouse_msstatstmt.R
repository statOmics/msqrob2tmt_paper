#### Mouse dataset: MSstatsTMT workflows ####

## This script models the mouse dataset using MSstatsTMT workflows.

####---- Setup environment ----####

# Read in the relevant libraries
library("MSstatsTMT")
library("BiocFileCache")

dataDir <- "data/"

source("scripts/utils.R") ## needed to get the hypotheses to test

####---- Load data ----####

## Download data (if needed)
bfc <- BiocFileCache()
psmFile <- bfcrpath(bfc, "https://zenodo.org/records/14767905/files/mouse_psms.txt?download=1")
annotFile <- bfcrpath(bfc, "https://zenodo.org/records/14767905/files/mouse_annotations.csv?download=1")

## Read data
mouse <- read.delim(psmFile)
annotation <- read.csv(annotFile)

####---- Run MSstatsTMT ----####

## Run MSstatsTMT on the full datasret
mouse <- PDtoMSstatsTMTFormat(
    input = mouse,
    annotation = annotation,
    which.proteinid = "Protein.Accessions",
    useUniquePeptide = TRUE, ## default TRUE
    rmPSM_withfewMea_withinRun = TRUE, ## default TRUE
    rmProtein_with1Feature = FALSE, ## default FALSE
    use_log_file = FALSE
)
mouse <- proteinSummarization(
    mouse, method = "msstats", # same as default
    global_norm = TRUE, # global peptide normalization between channels
    reference_norm = TRUE, # local protein normalization based on refernce channels
    MBimpute = TRUE,
    maxQuantileforCensored = NULL,
    remove_norm_channel = TRUE, # remove empty channels
    remove_empty_channel = TRUE, # remove norm channels
    use_log_file = FALSE
)

## Remove the Long_M samples, they are not used for modelling
mouse$ProteinLevelData <- mouse$ProteinLevelData[mouse$ProteinLevelData$Condition!= "Long_M", ]
mouse$ProteinLevelData$Condition <- as.character(mouse$ProteinLevelData$Condition)
mouse$FeatureLevelData <- mouse$FeatureLevelData[mouse$FeatureLevelData$Condition!= "Long_M", ]
mouse$FeatureLevelData$Condition <- as.character(mouse$FeatureLevelData$Condition)

####---- Statistical inference ----####

## hypothesesMouse is generated in `scripts/utils.R` to guarantee
## consistency with msqrob2TMT analyses
hypotheses <- hypothesesMouse
## Generate a contrast matrix from the hypotheses
contrasts <- t(msqrob2::makeContrast(
    hypotheses,
    parameterNames = c("ridgeConditionShort_HF", "ridgeConditionLong_HF", "ridgeConditionShort_LF", "ridgeConditionLong_LF")
))
dimnames(contrasts) <- lapply(dimnames(contrasts), function(x) gsub("ridgeCondition", "", x))

## Perform statistical inference (full dataset)
tests <- groupComparisonTMT(
    mouse,
    contrast.matrix = contrasts,
    moderated = TRUE
)$ComparisonResult

####---- Format and save output ----####

tests <- tests[, -8] ## remove the "issue" column
colnames(tests) <- c("Protein", "Comparison", "logFC", "se", "df", "pval", "adjPval")
tests$t <- NA ## for consistency with msqrob2TMT tables
tests$Model <- "MSstatsTMT"
tests$ModelVariable <- "Condition" ## for consistency with msqrob2TMT tables

## Store the modelled spikein1 data
saveRDS(data.frame(tests), paste0(dataDir, "mouse_model_MsstatsTMT.rds"))

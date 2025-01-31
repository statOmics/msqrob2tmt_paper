####---- MSstatsTMT workflow analyses on the UPS spike-in data ----####

## This script runs MSstatsTMT on the UPS spike-in data.

####---- Setup environment ----####

# Read in the relevant libraries
library("MSstatsTMT")
library("BiocFileCache")

dataDir <- "data/"

####---- Load data ----####

spikein1 <- readRDS(paste0(dataDir, "spikein1_input_msstatstmt.rds"))
## Retrieve the initial annotations
annotFileName <- bfcrpath(bfc, "https://zenodo.org/records/14767905/files/spikein1_annotations.csv?download=1")
annotation <- read.csv(annotFileName)

## Data set with no replication
spikein1NoRep <- spikein1[grepl("Mixture._01", spikein1$Spectrum.File), ]

####---- Run MSstatsTMT ----####

## Run MSstatsTMT on the full datasret
msstatsFull <- PDtoMSstatsTMTFormat(
    input = spikein1,
    annotation = annotation,
    which.proteinid = "Protein.Accessions",
    useUniquePeptide = TRUE, ## default TRUE
    rmPSM_withfewMea_withinRun = TRUE, ## default TRUE
    rmProtein_with1Feature = FALSE, ## default FALSE
    use_log_file = FALSE
) |>
    proteinSummarization(
        method = "msstats", # same as default
        global_norm = TRUE, # global peptide normalization between channels
        reference_norm = TRUE, # local protein normalization based on refernce channels
        MBimpute = TRUE,
        maxQuantileforCensored = NULL,
        remove_norm_channel = TRUE, # remove empty channels
        remove_empty_channel = TRUE, # remove norm channels
        use_log_file = FALSE
    )
## Repeat, but without technical replication
msstatsNoRep <- PDtoMSstatsTMTFormat(
    input = spikein1NoRep,
    annotation = annotation[annotation$Run %in% spikein1NoRep$Spectrum.File, ],
    which.proteinid = "Protein.Accessions",
    useUniquePeptide = TRUE, ## default TRUE
    rmPSM_withfewMea_withinRun = TRUE, ## default TRUE
    rmProtein_with1Feature = FALSE, ## default FALSE
    use_log_file = FALSE
) |>
    proteinSummarization(
        method = "msstats", # same as default
        global_norm = TRUE, # global peptide normalization between channels
        reference_norm = TRUE, # local protein normalization based on refernce channels
        MBimpute = TRUE,
        maxQuantileforCensored = NULL,
        remove_norm_channel = TRUE, # remove empty channels
        remove_empty_channel = TRUE, # remove norm channels
        use_log_file = FALSE
    )

####---- Statistical inference ----####

## Make contrast matrix
paramNames <- c("Condition1","Condition0.667","Condition0.5","Condition0.125")
contrasts <- combn(paramNames, 2, simplify = FALSE) |>
    sapply(function(x) paste(paste(x, collapse = " - "), " = 0")) |>
    msqrob2::makeContrast(paramNames) |>
    t()
dimnames(contrasts) <- lapply(dimnames(contrasts), function(x) gsub("Condition", "", x))

## Perform statistical inference (full dataset)
resultsFull <- groupComparisonTMT(
    msstatsFull,
    contrast.matrix = contrasts,
    moderated = TRUE
)$ComparisonResult
resultsFull$technicalReplicates <- TRUE

## Perform statistical inference (no technical replication)
resultsNoRep <- groupComparisonTMT(
    msstatsNoRep,
    contrast.matrix = contrasts,
    moderated = TRUE
)$ComparisonResult
resultsNoRep$technicalReplicates <- FALSE

out <- rbind(resultsFull, resultsNoRep)

####---- Format and save output ----####

out <- out[, -8] ## remove the "issue" column
colnames(out) <- c("Protein", "Comparison", "logFC", "se", "df", "pval", "adjPval", "technicalReplicates")
out$t <- NA
out$Model <- "MSstatsTMT"

## Store the modelled spikein1 data
saveRDS(data.frame(out), paste0(dataDir, "spikein1_model_MsstatsTMT.rds"))

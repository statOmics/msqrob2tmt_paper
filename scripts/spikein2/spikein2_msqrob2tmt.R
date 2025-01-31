
## msqrob2TMT workflow analyses on the spikein2 data ----####

## This script runs the msqrob2TMT workflows on the spikein2 data from
## O'Brien et al. 2024 where yeast proteins are spiked in a mouse
## background.
##
## We will perform msqrob2TMT analysis with the data preprocessed
## using O'Brien's code for custom preprocessing,as for all other
## methods, but we will also perform msqrob2TMT analysis on data
## further processed using msTrawler to decouple performance change
## associated with preprocessing and modelling.
##
## See also scripts/spikein2/spikein2_preprocess.R

####---- Setup environment ----####

library("msqrob2")
library("BiocFileCache")
library("tidyr")

## Load script with custom functions
source("scripts/utils.R")

dataDir <- "data/"

library("BiocParallel")
register(MulticoreParam(4)) ## limit analysis to 4 cores

## Sample annotations
sampleFile <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein2_annotations.csv?download=1") |>
    read.csv()
rownames(sampleFile) <- gsub("Sn", "Adjusted.Intensity", sampleFile$SampleID)
sampleFile$Run <- sub("_.*", "", sampleFile$SampleID)
sampleFile$Channel <- gsub(".*_|X|[.]Sn", "", sampleFile$SampleID)

####---- Get input data ----####

## These are the data used as input by all methods

## Read PSM data into a QFeatures object
spikein2 <- readRDS(paste0(dataDir, "spikein2_input.rds"))
spikein2 <- readQFeatures(
    spikein2,
    quantCols = grepl("Intensity",colnames(spikein2)),
    runCol = "Plex"
)
colData(spikein2) <- sampleFile
## Remove the normalisation channels
spikein2 <- subsetByColData(spikein2, spikein2$Bridge == 0)

## msqrob2TMT preprocessing
snames <- names(spikein2)
spikein2 <- logTransform(
    spikein2, snames, paste0(snames, "_log"), base = 2
)
spikein2 <- aggregateFeatures(
    spikein2,
    i = paste0(snames, "_log"),
    fcol = "Protein.ID",
    name = paste0(snames, "_proteins"),
    fun = MsCoreUtils::medianPolish,
    na.rm = TRUE
)
spikein2 <- joinAssays(
    spikein2, grep("log", names(spikein2)),
    name = "ions"
)
spikein2 <- joinAssays(
    spikein2, grep("proteins", names(spikein2)),
    name = "proteins"
)

####---- Get data preprocessed by msTrawler ----####

## This is the data that has been further processed by msTrawler (see
## scripts/spikein2/spikein2_preprocess.R).

## Get the data further preprocessed by msTrawler and convert them to
## a QFeatures object
spikein2MsTrawler <- paste0(dataDir, "spikein2_input_preprocessed.rds") |>
    readRDS() |>
    pivot_wider(
        id_cols = c("Protein.ID","Peptide","Plex"),
        names_from = Channel, values_from = lIntensity
    )
spikein2MsTrawler <- readQFeatures(
    spikein2MsTrawler,
    quantCols = grepl("X",colnames(spikein2MsTrawler)),
    runCol = "Plex"
)
sampleFile2 <- sampleFile
rownames(sampleFile2) <- gsub("Adjusted.Intensity", "Sn", rownames(sampleFile2))
colData(spikein2MsTrawler) <- sampleFile2
## Remove the normalisation channels
spikein2MsTrawler <- subsetByColData(
    spikein2MsTrawler, spikein2MsTrawler$Bridge == 0
)

## msqrob2TMT preprocessing
spikein2MsTrawler <- aggregateFeatures(
    spikein2MsTrawler,
    i = names(spikein2MsTrawler),
    fcol = "Protein.ID",
    name = paste0(names(spikein2MsTrawler), "_proteins"),
    fun = MsCoreUtils::medianPolish,
    na.rm = TRUE
)
spikein2MsTrawler <- joinAssays(
    spikein2MsTrawler, grep("\\d$", names(spikein2MsTrawler)),
    name = "ions"
)
spikein2MsTrawler <- joinAssays(
    spikein2MsTrawler, grep("proteins", names(spikein2MsTrawler)),
    name = "proteins"
)

####---- Define models and hypotheses ----####

## Define models
modelPsm <- ~ 0 + Dilution +
    Media +
    (1 | Run) +
    (1 | Run:Peptide) +
    (1 | Run:Channel)
modelProtein <- ~ 0 + Dilution +
    Media +
    (1 | Run)

## Create contrast matrix
paramNames <- unique(spikein2$Dilution) |>
    as.numeric() |>
    sort(decreasing = TRUE)
paramNames <- paste0("ridgeDilution", paramNames)
combinations <- sapply(
    combn(paramNames, 2, simplify = FALSE),
    function(x) paste(paste(x, collapse = " - "), " = 0")
)
L <- makeContrast(combinations, paramNames)


####---- Run msqrob2TMT workflows ----####

inference <- list()
for (preprocessing in c("msqrob2tmt", "msTrawler")) {
    if (preprocessing == "msqrob2tmt") {
        data <- spikein2
    } else {
        data <- spikein2MsTrawler
    }

    ####---- Modelling ----####

    ## msqrob2_rrilmm model
    data <- msqrob(
        data, i = "proteins",
        formula = modelProtein,
        modelColumnName = "msqrob2_rrilmm",
        ridge = TRUE, robust = TRUE
    )

    ## msqrob2_psm_rrilmm model
    data <- msqrobAggregate( ## modelling
        data, i = "ions",
        formula = modelPsm,
        fcol = "Protein.ID",
        modelColumnName = "msqrob2_psm_rrilmm",
        name = "proteins_msqrob",
        ridge = TRUE, robust = TRUE,
        aggregateFun = colSums
    )
    ## msqrob2_psm_rrilmm_refit model
    ## Add a new assay to avoid overwriting the msqrob2_psm_rrilmm results
    data <- addAssay(
        data, data[["proteins_msqrob"]], "proteins_msqrob_refit"
    )
    ## Retrieve one-hit wonders
    counts <- aggcounts(data[["proteins_msqrob"]])
    oneHitProteins <- rownames(counts)[rowMax(counts) == 1]
    ## Remove the channel within run random effect
    simplifiedModel <- as.character(modelPsm) |>
        paste(collapse = " ") |>
        sub(pattern = " \\+ \\(1 \\| Run:Channel\\)", replacement = "") |>
        as.formula()
    ## Refit
    data <- msqrobRefit(
        data, i = "ions",
        subset = oneHitProteins,
        fcol = "Protein.ID",
        formula = simplifiedModel,
        ridge = TRUE,
        robust = TRUE,
        name = "proteins_msqrob_refit",
        modelColumnName = "msqrob2_psm_rrilmm",
        aggregateFun = colSums
    )

    ####---- Statistical inference ----####

    data <- hypothesisTest(
        data, i = "proteins", L,
        modelColumn = "msqrob2_rrilmm"
    )
    data <- hypothesisTest( ## Statistical inference
        data, i = "proteins_msqrob", L,
        modelColumn = "msqrob2_psm_rrilmm"
    )
    data <- hypothesisTest(
        data, i = "proteins_msqrob_refit", contrast = L,
        modelColumn = "msqrob2_psm_rrilmm"
    )

    ####---- Retrieve results ----####

    out <- c("proteins", "proteins_msqrob", "proteins_msqrob_refit") |>
        lapply(function(i) {
            rd <- rowData(data)[[i]]
            tests <- rd[, colnames(L)]
            tests <- lapply(names(tests), function(ii) {
                test <- tests[[ii]]
                test$Comparison <- sub(".*(ridge)?Dilution(.*) - (ridge)?Dilution(.*)", "\\2 - \\4", ii)
                test$Protein <- rownames(test)
                test
            })
            tests <- do.call(rbind, tests)
            modelName <- colnames(rd)[grepl("msqrob", colnames(rd))]
            modelName <- ifelse(grepl("refit", i), paste0(modelName, "_refit"), modelName)
            tests$Model <- modelName
            tests
        }) |>
        do.call(what = rbind)
    out$Preprocessing <- preprocessing
    inference[[preprocessing]] <- out
}

####---- Save output ----####

inference <- do.call(rbind, inference)
saveRDS(inference, paste0(dataDir, "spikein2_model_msqrob2tmt.rds"))

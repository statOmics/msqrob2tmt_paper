####---- Compare preprocessing workflows ----####

## This script runs MSstatsTMT and msqrob2TMT workflows on the UPS
## spike-in data preprocessed using 4 strategies:
## 1. imp-refnorm: perform imputation and reference normalisation
## 2. noimp-refnorm: perform no imputation and reference normalisation.
## 3. imp-refnorm: perform imputation and no reference normalisation.
## 3. noimp-refnorm: perform no imputation and no reference normalisation.

####---- Setup environment ----####

## Read in the relevant libraries
library("msqrob2")
library("QFeatures")
library("MSstatsTMT")
library("tidyr")

dataDir <- "data/"

library("BiocParallel")
register(MulticoreParam(4))

####---- Load data ----####

spikein1 <- readRDS(paste0(dataDir, "spikein1_input_msstatstmt.rds"))
## Retrieve the initial annotations
library("BiocFileCache")
annotation <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein1_annotations.csv?download=1") |>
    read.csv()
annotation$Run <- gsub("[.]", "", annotation$Run)
annotation$runCol <- annotation$Run
annotation$quantCols <- annotation$Channel
## Data set with no replication
spikein1 <- spikein1[grepl("Mixture._01", spikein1$Spectrum.File), ]

## Make MSstatsTMT required format
spikein1 <- PDtoMSstatsTMTFormat(
    input = spikein1,
    annotation = annotation,
    which.proteinid = "Protein.Accessions", ## same as default
    use_log_file = FALSE
)

####---- Make contrasts ----####

paramNames <- c("Condition1","Condition0.667","Condition0.5","Condition0.125")

## Make contrast matrix for MSstatsTMT
contrastsMSstats <- combn(paramNames, 2, simplify = FALSE) |>
    sapply(function(x) paste(paste(x, collapse = " - "), " = 0")) |>
    msqrob2::makeContrast(paramNames) |>
    t()
dimnames(contrastsMSstats) <- lapply(dimnames(contrastsMSstats), function(x) gsub("Condition", "", x))

## Make contrast matrix for MSstatsTMT
combinations <- sapply(
    combn(paramNames, 2, simplify = FALSE),
    function(x) paste(paste(x, collapse = " - "), " = 0")
)
contrastsMsqrob <- makeContrast(combinations, paramNames)


####---- Run models ----####

strategies <- expand.grid(
    imputation = c("imp", "noimp"),
    normalisation = c("refnorm", "norefnorm")
)

out <- list()
for (i in 1:nrow(strategies)) {
    imputation <- strategies$imputation[i]
    normalisation <- strategies$normalisation[i]

    cat("Running the analysis with params:", imputation, "and", normalisation, "\n")

    MBimpute <- ifelse(imputation == "imp", TRUE, FALSE)
    if (normalisation == "norefnorm") {
        data <- spikein1[spikein1$Condition != "Norm",]
        data$Condition <- factor(data$Condition)
        reference_norm <- FALSE
    } else {
        data <- spikein1
        reference_norm <- TRUE
    }

    ## Data preprocessing
    data <- proteinSummarization(
        data,
        method = "msstats", # same as default
        global_norm = TRUE, # global peptide normalization between channels
        reference_norm = reference_norm, # local protein normalization based on refernce channels
        MBimpute = MBimpute,
        maxQuantileforCensored = NULL,
        remove_norm_channel = TRUE, # remove empty channels
        remove_empty_channel = TRUE, # remove norm channels
        use_log_file = FALSE
    )

    ## Run MSstatsTMT inference
    resultMsStatsTMT <- groupComparisonTMT(
        data, contrast.matrix = contrastsMSstats, moderated = TRUE
    )$ComparisonResult
    resultMsStatsTMT <- resultMsStatsTMT[, -8] ## remove the "issue" column
    colnames(resultMsStatsTMT) <- c("Protein", "Comparison", "logFC", "se", "df", "pval", "adjPval")
    resultMsStatsTMT$t <- NA
    resultMsStatsTMT$Model <- paste0("MSstatsTMT_", imputation, "_", normalisation)

    ## Format to QFeatures
    data <- pivot_wider(
        data$ProteinLevelData,
        id_cols = c("Mixture", "TechRepMixture", "Run", "Protein"),
        names_from = "Channel",
        values_from ="Abundance"
    )
    data  <- readQFeatures(
        data,
        colData = annotation,
        runCol = "Run",
        quantCols = grep("^1", colnames(data))
    )
    for (ii in names(data)) {
        rownames(data[[ii]]) <- rowData(data[[ii]])$Protein
    }
    data <- joinAssays(data, names(data), "proteins")

    ## Run msqrob2TMT
    data <- msqrob(
        data, i = "proteins",
        formula = ~ 0 + Condition + (1 | runCol),
        ridge = FALSE,
        robust = FALSE
        )
    data <- hypothesisTest(data, i = "proteins", contrast = contrastsMsqrob)

    ## Retrieve and format msqrob2TMT results
    tests <- rowData(data[["proteins"]])[, colnames(contrastsMsqrob)]
    tests <- lapply(names(tests), function(ii) {
        test <- tests[[ii]]
        test$Comparison <- sub(".*(ridge)?Condition(.*) - (ridge)?Condition(.*)", "\\2 - \\4", ii)
        test$Protein <- rownames(test)
        test
    })
    resultMssqrob2TMT <- do.call(rbind, tests)
    resultMssqrob2TMT$Model <- paste0("msqrob2tmt_lmm_", imputation, "_", normalisation)

    out <- c(out, list(resultMsStatsTMT, resultMssqrob2TMT))
}

out <- do.call(rbind, out)
out$technicalReplicates <- FALSE
saveRDS(out, paste0(dataDir, "spikein1_model_compare_preprocessing.rds"))

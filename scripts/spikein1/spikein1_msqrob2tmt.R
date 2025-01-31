
## msqrob2TMT workflow analyses on the UPS spike-in data ----####

## This script runs the msqrob2TMT workflows on the UPS spike-in data.

####---- Setup environment ----####

library("msqrob2")

## Load script with custom functions
source("scripts/utils.R")

dataDir <- "data/"

library("BiocParallel")
register(MulticoreParam(4)) ## limit analysis to 4 cores

####---- Load data ----####

spikein1 <- readRDS(paste0(dataDir, "spikein1_input_msqrob2tmt.rds"))
## Remove reference channel (not used by msqrob2)
spikein1 <- subsetByColData(spikein1, spikein1$Condition != "Norm")

####---- Preprocessing ----####

sNames <- names(spikein1)[1:15]
## Normalisation
spikein1 <- normalize(
    spikein1, paste0(sNames, "_log"), name = paste0(sNames, "_norm"),
    method = "center.median"
)
## Summarisation
spikein1 <- aggregateFeatures(
    spikein1, i = paste0(sNames, "_norm"), name = paste0(sNames, "_proteins"),
    fcol = "Protein.Accessions", fun = MsCoreUtils::medianPolish,
    na.rm=TRUE
)
## Join sets
spikein1 <- joinAssays(spikein1, paste0(sNames, "_log"), "ions")
spikein1 <- joinAssays(spikein1, paste0(sNames, "_norm"), "ions_norm")
spikein1 <- joinAssays(spikein1, paste0(sNames, "_proteins"), "proteins")

####---- Define models ----####

## Generate a table with all arguments needed to run the different
## flavours of msqrob2TMT workflows
## lmm = LMM with no ridge, no robust
## rlmm = LMM with no ridge, robust
## rilmm = LMM with ridge, no robust
## rrilmm = LMM with ridge and robust

## Protein-level models
proteinModelArgs <- data.frame(
    model = c("msqrob2_lmm", "msqrob2_rlmm", "msqrob2_rilmm", "msqrob2_rrilmm"),
    robust = c(FALSE, TRUE, FALSE, TRUE),
    ridge = c(FALSE, FALSE, TRUE, TRUE)
)
proteinModelArgs$formula <- "~ 0 + Condition + (1|Run)"
proteinModelArgs$i <- "proteins"
proteinModelArgs$psmLevel <- FALSE
## PSM-level models
psmModelArgs <- proteinModelArgs
psmModelArgs$model <- sub("msqrob2_", "msqrob2_psm_", psmModelArgs$model)
psmModelArgs$formula <- paste0(psmModelArgs$formula, " + (1 | Run:ionID) + (1 | Run:Channel)")
psmModelArgs$i <- "ions_norm"
psmModelArgs$psmLevel <- TRUE
## Combine
modelArgs <- rbind(psmModelArgs, proteinModelArgs)
## Duplicate these parameters to run the models on the full data set
## or without technical replicates
modelArgs$technicalReplicates <- FALSE
modelArgsFull <- modelArgs[c(4, 8), ]
modelArgsFull$technicalReplicates <- TRUE
modelArgsFull$formula <- paste0(modelArgsFull$formula, " + (1 | Mixture)")

## Combine all model arguments
modelArgs <- rbind(modelArgs, modelArgsFull)

####---- Run models and statistical inference ----####

## Use the table to run the different msqrob2TMT workflows
out <- lapply(1:nrow(modelArgs), function(ii) {
    args <- modelArgs[ii, ]
    cat("Running model:", args$model, "with",
        ifelse(args$technicalReplicates, "", "no"), "replication.\n")
    do.call(msqrob2Workflow, c(data = spikein1, args))
})

## Refit the msqrob2_psm_rrilmm models to assess impact of fitting
## one-hit-wonder proteins (done for both the full data set and the
## dataset without technical replication)
for (ii in grep("psm_rrilmm", modelArgs$model)) {
    args <- modelArgs[ii, ]
    data <- do.call(runMsqrob2, c(data = spikein1, args))
    ## Retrieve one-hit wonders
    counts <- aggcounts(data[["proteins_msqrob2"]])
    oneHitProteins <- rownames(counts)[rowMax(counts) == 1]
    ## Refit happens here
    simplifiedModel <- sub(" \\+ \\(1 \\| Run:Channel\\)", "", args$formula)
    data <- msqrobRefit(
        data,
        i = args$i,
        subset = oneHitProteins,
        fcol = "Protein.Accessions",
        formula = as.formula(simplifiedModel),
        ridge = TRUE,
        robust = TRUE,
        name = "proteins_msqrob2",
        modelColumnName = "msqrobModels",
        aggregateFun = colSums
    )
    ## Run inference
    L <- generateContrast(args$ridge)
    data <- hypothesisTest(
        data, i = "proteins_msqrob2", contrast = L,
        resultsColumnNamePrefix = "test_"
    )
    ## Retrieve results
    refitResults <- inferenceResults(data, i = "proteins_msqrob2")
    refitResults$Model <- paste0(args$model, "_refit")
    refitResults$technicalReplicates <- args$technicalReplicates
    out <- c(out, list(refitResults))
}

####---- Save output ----####

## Combine all results in a single table
out <- do.call(rbind, out)
## Save
saveRDS(out, paste0(dataDir, "spikein1_model_msqrob2tmt.rds"))

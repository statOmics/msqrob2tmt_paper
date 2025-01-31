
#### Mouse dataset: msqrob2TMT workflows ####

## This script models the mouse dataset using msqrob2TMT workflows,
## modelling both at the PSM and at the protein level. For the
## protein-level modelling, data are summarised for each fraction
## separately. See the `mouse_msqrob2tmt_mixture.R` script to
## summarise proteins by combining peptides across the different
## fractions within a mixture.

####---- Setup environment ----####

## Read in the relevant libraries
library("msqrob2")
library("dplyr")
library("BiocFileCache")

library("BiocParallel")
register(MulticoreParam(4)) ## limit analysis to 4 cores

dataDir <- "data/"

source("scripts/utils.R") ## needed for msqrobRefit()

####---- Load the data ----####

## Download data (if needed)
bfc <- BiocFileCache()
psmFile <- bfcrpath(bfc, "https://zenodo.org/records/14767905/files/mouse_psms.txt?download=1")
annotFile <- bfcrpath(bfc, "https://zenodo.org/records/14767905/files/mouse_annotations.csv?download=1")

## Read data
psms <- read.delim(psmFile)
coldata <- read.csv(annotFile)
coldata$File.Name <- coldata$Run
coldata$Duration <- gsub("_.*", "", coldata$Condition)
coldata$Diet <- gsub(".*_", "", coldata$Condition)

## Format as a QFeatures
coldata$runCol <- coldata$Run
coldata$quantCols <- paste0("Abundance..", coldata$Channel)
mouse <- readQFeatures(psms, colData = coldata,
                       quantCols = unique(coldata$quantCols),
                       runCol = "Spectrum.File", name = "psms")

####---- PSM filtering ----####

## Remove failed or ambiguous protein inference
mouse <- filterFeatures(
    mouse, ~ Protein.Accessions != "" & ## Remove failed protein inference
        !grepl(";", Protein.Accessions)) ## Remove protein groups

## Missing data filtering
mouse <- zeroIsNA(mouse, names(mouse))
mouse <- filterNA(mouse, names(mouse), pNA = 0.7)

## Remove PSMs that map to multiple
for (i in names(mouse)) {
    rowdata <- rowData(mouse[[i]])
    rowdata$ionID <- paste0(rowdata$Annotated.Sequence, rowdata$Charge)
    rowdata$rowsums <- rowSums(assay(mouse[[i]]), na.rm=TRUE)
    rowdata <- data.frame(rowdata) |>
        group_by(ionID) |>
        mutate(psmRank = rank(-rowsums))
    rowData(mouse[[i]]) <- DataFrame(rowdata)
}
mouse <- filterFeatures(mouse, ~ psmRank == 1)
mouse <- filterFeatures(mouse, ~ rowsums > 0)

## Redefine the rownames of each set using on the ion identifier as it
## is now unique within each run.
for(i in names(mouse)){
    rownames(mouse[[i]]) <- rowData(mouse[[i]])$ionID
}

####----  Preprocessing ----####

sNames <- names(mouse)
mouse <- logTransform(
    mouse, sNames, name = paste0(sNames, "_log"), base = 2
)
mouse <- normalize(
    mouse, paste0(sNames, "_log"), name = paste0(sNames, "_norm"),
    method = "center.median"
)
mouse <- aggregateFeatures(
    mouse, i = paste0(sNames, "_norm"), name = paste0(sNames, "_proteins"),
    fcol = "Protein.Accessions", fun = MsCoreUtils::medianPolish,
    na.rm=TRUE
)
mouse <- subsetByColData(
    mouse, mouse$Condition != "Norm" & mouse$Condition != "Long_M"
)
mouse <- joinAssays(mouse, paste0(sNames, "_norm"), "ions_norm")
mouse <- joinAssays(mouse, paste0(sNames, "_proteins"), "proteins")

####---- Simulate mock conditions ----####

mouseMock <- mouse
## Remove a mixture to avoid imbalanced data issues
mouseMock <- subsetByColData(
    mouseMock,
    mouseMock$Condition != "Long_HF" &
        mouseMock$Mixture != "PAMI-194_Mouse_U-Dd"
)
## Assign the mock conditions
mouseMock$Mock[mouseMock$Mixture == "PAMI-176_Mouse_A-J"] <- c("A","B","A","B","A","B")
mouseMock$Mock[mouseMock$Mixture == "PAMI-176_Mouse_K-T"] <- c("A","B","B","A","B","A")

## Check the randomisation
table(Condition = mouseMock$Condition, Mock = mouseMock$Mock)
table(Mixture = mouseMock$Mixture, Mock = mouseMock$Mock)

####---- Define models ----####

## Define models
modelProtein <- ~ 0 +  ## Remove intercept to simplify parameter intrepretation
    Condition + ## fixed effect for experimental condition
    (1 | Mixture) + ## random effect for mixture
    (1 | Run) + ## random effect for MS run
    (1 | BioReplicate)  ## random effect for biological replicate (mouse)
modelPsm <- ~ 0 +  ## Remove intercept to simplify parameter intrepretation
    Condition + ## fixed effect for experimental condition
    (1 | Mixture) + ## random effect for mixture
    (1 | Run) + ## random effect for MS run
    (1 | Run:ionID) + ## random effect for spectrum, i.e. ionID nested in run
    (1 | Run:Channel) + ## random effect for channel nested in MS run
    (1 | BioReplicate)  ## random effect for biological replicate (mouse)


## hypothesesMouse is generated in `scripts/utils.R` to guarantee
## consistency with MSstatsTMT analyses
hypotheses <- hypothesesMouse

## Generate a contrast matrix from the hypotheses
L <- makeContrast(
    hypotheses,
    parameterNames = c("ridgeConditionShort_HF", "ridgeConditionLong_HF", "ridgeConditionShort_LF", "ridgeConditionLong_LF")
)

## Same but with the Mock condition
modelProteinMock <- ~ 0 +  ## Remove intercept to simplify parameter intrepretation
    Mock + ## fixed effect for the mock condition
    Condition + ## fixed effect for experimental condition
    (1 | Mixture) + ## random effect for mixture
    (1 | Run) + ## random effect for MS run
    (1 | BioReplicate)  ## random effect for biological replicate (mouse)
modelPsmMock <- ~ 0 +  ## Remove intercept to simplify parameter intrepretation
    Mock + ## fixed effect for the mock condition
    Condition + ## fixed effect for experimental condition
    (1 | Mixture) + ## random effect for mixture
    (1 | Run) + ## random effect for MS run
    (1 | Run:ionID) + ## random effect for spectrum, i.e. ionID nested in run
    (1 | Run:Channel) + ## random effect for channel nested in MS run
    (1 | BioReplicate)  ## random effect for biological replicate (mouse)
## Generate a contrast matrix for mock condition
LMock <- makeContrast(
    "ridgeMockA - ridgeMockB = 0",
    parameterNames = c("ridgeMockA", "ridgeMockB")
)

####---- Run msqrob2 for real condition ----####

## msqrob2_rrilmm model
mouse <- msqrob(
    mouse, i = "proteins",
    formula = modelProtein,
    ridge = TRUE, robust = TRUE
)
## Statistical inference
mouse <- hypothesisTest(mouse, i = "proteins", L)

## msqrob2_psm_rrilmm model
mouse <- msqrobAggregate( ## modelling
    mouse, i = "ions_norm",
    formula = modelPsm,
    fcol = "Protein.Accessions",
    name = "proteins_msqrob",
    ridge = TRUE, robust = TRUE,
    aggregateFun = colMedians
)
## Statistical inference
mouse <- hypothesisTest(mouse, i = "proteins_msqrob", L)

## msqrob2_psm_rrilmm_refit model
## Add a new assay to avoid overwriting the msqrob2_psm_rrilmm results
mouse <- addAssay(
    mouse, mouse[["proteins_msqrob"]], "proteins_msqrob_refit"
)
## Retrieve one-hit wonders
counts <- aggcounts(mouse[["proteins_msqrob_refit"]])
oneHitProteins <- rownames(counts)[rowMax(counts) == 1]
## Remove the channel wihtin run random effect
simplifiedModel <- as.character(modelPsm) |>
    paste(collapse = " ") |>
    sub(pattern = " \\+ \\(1 \\| Run:Channel\\)", replacement = "") |>
    as.formula()
## Refit
mouse <- msqrobRefit(
    mouse, i = "ions_norm",
    subset = oneHitProteins,
    fcol = "Protein.Accessions",
    formula = simplifiedModel,
    ridge = TRUE,
    robust = TRUE,
    name = "proteins_msqrob_refit",
    modelColumnName = "msqrobModels",
    aggregateFun = colSums
)
mouse <- hypothesisTest(
    mouse, i = "proteins_msqrob_refit", L, overwrite = TRUE
)

####---- Run msqrob2 for mock condition ----####

## msqrob2_rrilmm model
mouseMock <- msqrob(
    mouseMock, i = "proteins",
    formula = modelProteinMock,
    ridge = TRUE, robust = TRUE
)
mouseMock <- hypothesisTest( ## Statistical inference
    mouseMock, i = "proteins", LMock
)

## msqrob2_psm_rrilmm model
mouseMock <- msqrobAggregate( ## modelling
    mouseMock, i = "ions_norm",
    formula = modelPsmMock,
    fcol = "Protein.Accessions",
    name = "proteins_msqrob",
    ridge = TRUE, robust = TRUE,
    aggregateFun = colMedians
)
mouseMock <- hypothesisTest( ## Statistical inference
    mouseMock, i = "proteins_msqrob", LMock
)

## msqrob2_psm_rrilmm_refit model
## Add a new assay to avoid overwriting the msqrob2_psm_rrilmm results
mouseMock <- addAssay(
    mouseMock, mouseMock[["proteins_msqrob"]], "proteins_msqrob_refit"
)
## Retrieve one-hit wonders
counts <- aggcounts(mouseMock[["proteins_msqrob_refit"]])
oneHitProteins <- rownames(counts)[rowMax(counts) == 1]
## Remove the channel wihtin run random effect
simplifiedModel <- as.character(modelPsm) |>
    paste(collapse = " ") |>
    sub(pattern = " \\+ \\(1 \\| Run:Channel\\)", replacement = "") |>
    as.formula()
## Refit
mouseMock <- msqrobRefit(
    mouseMock, i = "ions_norm",
    subset = oneHitProteins,
    fcol = "Protein.Accessions",
    formula = simplifiedModel,
    ridge = TRUE,
    robust = TRUE,
    modelColumnName = "msqrobModels",
    name = "proteins_msqrob_refit",
    aggregateFun = colSums
)
mouseMock <- hypothesisTest(
    mouseMock, i = "proteins_msqrob_refit", LMock, overwrite = TRUE
)

####---- Retrieve the results ----####

modelNames <- c(
    proteins = "msqrob2_rrilmm",
    proteins_msqrob = "msqrob2_psm_rrilmm",
    proteins_msqrob_refit = "msqrob2_psm_rrilmm_refit"
)

out <- lapply(c("Condition", "Mock"), function(j) {
    data <- if (j == "Condition") mouse else mouseMock
    testNames <- if (j == "Condition") colnames(L) else colnames(LMock)
    pattern <- paste0("(ridge)?", j)
    outj <- lapply(1:length(modelNames), function(i) {
        rd <- rowData(data[[names(modelNames)[i]]])
        tests <- rd[, testNames, drop = FALSE]
        ## A little cleaning
        tests <- lapply(names(tests), function(ii) {
            test <- tests[[ii]]
            test$Comparison <- gsub(pattern, "", ii)
            test$Protein <- rownames(test)
            test
        })
        ## Return a single table for all comparisons
        tests <- do.call(rbind, tests)
        tests$Model <- modelNames[[i]]
        tests
    })
    outj <- do.call(rbind, outj)
    outj$ModelVariable <- j
    outj
})

####---- Save output ----####

## Combine all results in a single table
out <- do.call(rbind, out)
## Save
saveRDS(out, paste0(dataDir, "mouse_model_msqrob2tmt.rds"))

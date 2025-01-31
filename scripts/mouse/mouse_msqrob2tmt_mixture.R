
#### Mouse dataset: msqrob2TMT protein workflow by summarising mixtures ####

## This script models the mouseMix dataset using msqrob2TMT workflows,
## modelling both at the protein level. The summarisation will combine
## peptides found across the multiple fractions (MS runs) within a
## mixture. Note that this approach is the same approach as
## implemented in MSstatsTMT. See the `mouseMix_msqrob2tmt.R` script to
## summarise the peptides for each fraction separately.

####---- Setup environment ----####

## Read in the relevant libraries
library("msqrob2")
library("dplyr")
library("BiocFileCache")

library("BiocParallel")
register(MulticoreParam(4)) ## limit analysis to 4 cores

dataDir <- "data/"

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
coldata <- coldata |>
    select(-Run,-Fraction,-File.Name) |>
    unique() |>
    mutate(runCol = Mixture,
           quantCols = paste0("Abundance..", Channel)
    )
psms$Mixture <- sub("(.*)_TMT.*", "\\1", psms$Spectrum.File)
mouseMix <- readQFeatures(psms, colData = coldata,
                          quantCols = unique(coldata$quantCols),
                          runCol = "Mixture")

####---- PSM filtering ----####

## Remove failed or ambiguous protein inference
mouseMix <- filterFeatures(
    mouseMix, ~ Protein.Accessions != "" & ## Remove failed protein inference
        !grepl(";", Protein.Accessions)) ## Remove protein groups

## Missing data filtering
mouseMix <- zeroIsNA(mouseMix, names(mouseMix))
mouseMix <- filterNA(mouseMix, names(mouseMix), pNA = 0.7)

## Remove PSMs that map to multiple
for (i in names(mouseMix)) {
    rowdata <- rowData(mouseMix[[i]])
    rowdata$ionID <- paste0(rowdata$Annotated.Sequence, rowdata$Charge)
    rowdata$rowsums <- rowSums(assay(mouseMix[[i]]), na.rm=TRUE)
    rowdata <- data.frame(rowdata) |>
        group_by(ionID) |>
        mutate(psmRank = rank(-rowsums))
    rowData(mouseMix[[i]]) <- DataFrame(rowdata)
}
mouseMix <- filterFeatures(mouseMix, ~ psmRank == 1)
mouseMix <- filterFeatures(mouseMix, ~ rowsums > 0)

## Redefine the rownames of each set using on the ion identifier as it
## is now unique within each run.
for(i in names(mouseMix)){
    rownames(mouseMix[[i]]) <- rowData(mouseMix[[i]])$ionID
}

####----  Preprocessing ----####

sNames <- names(mouseMix)
mouseMix <- logTransform(
    mouseMix, sNames, name = paste0(sNames, "_log"), base = 2
)
mouseMix <- normalize(
    mouseMix, paste0(sNames, "_log"), name = paste0(sNames, "_norm"),
    method = "center.median"
)
mouseMix <- aggregateFeatures(
    mouseMix, i = paste0(sNames, "_norm"), name = paste0(sNames, "_proteins"),
    fcol = "Protein.Accessions", fun = MsCoreUtils::medianPolish,
    na.rm=TRUE
)
mouseMix <- subsetByColData(
    mouseMix, mouseMix$Condition != "Norm" & mouseMix$Condition != "Long_M"
)
mouseMix <- joinAssays(mouseMix, paste0(sNames, "_proteins"), "proteins")

####---- Simulate mock conditions ----####

mouseMixMock <- mouseMix
## Remove a mixture to avoid imbalanced data issues
mouseMixMock <- subsetByColData(
    mouseMixMock,
    mouseMixMock$Condition != "Long_HF" &
        mouseMixMock$Mixture != "PAMI-194_Mouse_U-Dd"
)
## Assign the mock conditions
mouseMixMock$Mock[mouseMixMock$Mixture == "PAMI-176_Mouse_A-J"] <- c("A","B","A","B","A","B")
mouseMixMock$Mock[mouseMixMock$Mixture == "PAMI-176_Mouse_K-T"] <- c("A","B","B","A","B","A")
## Check the randomisation
table(Condition = mouseMixMock$Condition, Mock = mouseMixMock$Mock)
table(Mixture = mouseMixMock$Mixture, Mock = mouseMixMock$Mock)

####---- Define models ----####

## Define models
modelMixture <- ~ 0 + ## Remove intercept to simplify parameter intrepretation
    Condition + ## fixed effect for Condition with interaction
    (1 | Mixture)  ## random effect for mixture

## hypothesesMouse is generated in `scripts/utils.R` to guarantee
## consistency with MSstatsTMT analyses
hypotheses <- hypothesesMouse

## Generate a contrast matrix from the hypotheses
L <- makeContrast(
    hypotheses,
    parameterNames = c("ridgeConditionShort_HF", "ridgeConditionLong_HF", "ridgeConditionShort_LF", "ridgeConditionLong_LF")
)

## Same but for the Mock condition
modelMixtureMock <- ~ 0 + ## Remove intercept to simplify parameter intrepretation
    Mock + ## fixed effect for the mock condition
    Condition + ## fixed effect for experimental condition
    (1 | Mixture)
## Generate a contrast matrix for mock condition
LMock <- makeContrast(
    "ridgeMockA - ridgeMockB = 0",
    parameterNames = c("ridgeMockA", "ridgeMockB")
)

####---- Run msqrob2 for real condition ----####

## msqrob2_rrilmm model
mouseMix <- msqrob(
    mouseMix, i = "proteins",
    formula = modelMixture,
    ridge = TRUE, robust = TRUE
)
mouseMix <- hypothesisTest( ## Statistical inference
    mouseMix, i = "proteins", L
)

####---- Run msqrob2 for mock condition ----####

## msqrob2_rrilmm model
mouseMixMock <- msqrob(
    mouseMixMock, i = "proteins",
    formula = modelMixtureMock,
    ridge = TRUE, robust = TRUE
)
mouseMixMock <- hypothesisTest( ## Statistical inference
    mouseMixMock, i = "proteins", LMock
)

####---- Retrieve the results ----####

out <- lapply(c("Condition", "Mock"), function(j) {
    data <- if (j == "Condition") mouseMix else mouseMixMock
    testNames <- if (j == "Condition") colnames(L) else colnames(LMock)
    pattern <- paste0("(ridge)?", j)
    rd <- rowData(data[["proteins"]])
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
    tests$Model <- "msqrob2_rrilmm_mixture"
    tests$ModelVariable <- j
    tests
})

####---- Save output ----####

## Combine all results in a single table
out <- do.call(rbind, out)
## Save
saveRDS(out, paste0(dataDir, "mouse_model_msqrob2tmt_mixture.rds"))


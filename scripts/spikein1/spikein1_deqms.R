####---- DEqMS workflow analyses on the UPS spike-in data ----####

## This script runs DEqMS on the UPS spike-in data.

####---- Setup environment ----####

# Read in the relevant libraries
library("DEqMS")
library("BiocFileCache")

dataDir <- "data/"

####---- Load data ----####

## Read PSM data
spikein1 <- readRDS(paste0(dataDir, "spikein1_input_deqms.rds"))
## Read sample annotations
annotations <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein1_annotations.csv?download=1") |>
    read.csv()
## Match PSM data and annotations
row.names(annotations) <- paste0(annotations$Run, "_Abundance..", annotations$Channel)
annotations <- annotations[colnames(spikein1)[-(1:2)], ]

## Get PSM counts per protein for each run
counts <- lapply(split(colnames(spikein1)[-(1:2)], annotations$Run), function(ii) {
    observed <- rowSums(!is.na(spikein1[, ii])) != 0
    proteinsInRun <- spikein1$Protein.Accessions[observed]
    data.frame(table(Protein.Accessions = proteinsInRun))
})
counts <- Reduce(
    function(x, y) merge(x, y, by = "Protein.Accessions", all = TRUE),
    counts
)
countsn <- counts$Protein.Accessions
counts <- rowMins(as.matrix(counts[, -1]), na.rm = TRUE)
names(counts) <- countsn

####---- DEqMS preprocessing ----####

## medianSweeping() = median sweep summarisation + col median
## normalisation
spikein1 <- medianSweeping(spikein1, group_col = 2)

####---- Run DEqMS ----####

## Make design matrix
## Fixed effect for condition and fixed effect for run
design <- model.matrix(~ 0 + Condition + Run, data = annotations)

## Make contrasts
contrast <- colnames(design) |>
    grep(pattern = "Condition", value = TRUE) |>
    sort(decreasing = TRUE) |>
    combn(2, paste, collapse = " - ")
contrast <- makeContrasts(contrasts = contrast, levels = design)

## Run DEqMS modelling
fit1 <- lmFit(spikein1, design)
## Note that some proteins have very low df, we will remove them as
## this will lead to data inconsistencies downstream the script
fit1 <- lmFit(spikein1[fit1$df.residual > 1, ], design)
fit2 <- eBayes(contrasts.fit(fit1, contrasts = contrast))
fit2$count <- counts[rownames(fit2$coefficients)]
fit3 <- spectraCounteBayes(fit2)

####---- Format and save output ----####

out <- lapply(colnames(fit3$coefficients), function(ii) {
    test <- outputResult(fit3, coef_col = ii)
    test <- test[, c("logFC","t","sca.adj.pval","sca.P.Value")]
    colnames(test) <- c("logFC","t","adjPval","pval")
    test$Protein <- rownames(test)
    test$df <- NA
    test$se <- NA
    test$Comparison <- gsub("Condition", "", ii)
    test$Model <- "DEqMS"
    test$technicalReplicates <- FALSE
    test
})
out <- do.call(rbind, out)

## Store the modelled spikein1 data
saveRDS(out, paste0(dataDir, "spikein1_model_DEqMS.rds"))

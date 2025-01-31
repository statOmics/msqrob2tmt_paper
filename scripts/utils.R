
####---- Hypotheses for Mouse data ----####

## We define the hypotheses to test for the mouse data set:
## 1. Difference between low fat and high fat diet after short
## duration
## 2. Difference between low fat and high fat diet after long duration
## 3. Average difference between low fat and high fat diet
## 4. Interaction: does the diet effect change according to duration?
## Note that storing these hypotheses in this script ensures that both
## msqrob2TMT and MSstatsTMT workflows perform the same hypothesis
## testing
hypothesesMouse <- c(
    "ridgeConditionShort_LF - ridgeConditionShort_HF = 0", ## 1.
    "ridgeConditionLong_LF - ridgeConditionLong_HF = 0", ## 2.
    "(ridgeConditionShort_LF + ridgeConditionLong_LF)/2 - (ridgeConditionShort_HF + ridgeConditionLong_HF)/2 = 0", ## 3.
    "(ridgeConditionShort_LF - ridgeConditionShort_HF) - (ridgeConditionLong_LF - ridgeConditionLong_HF) = 0" ## 4.
)


####---- msqrob2 automation utils ----####

msqrobRefit <- function(object, formula, i, subset, fcol, name,
                        modelColumnName, ...) {
    seti <- getWithColData(object, i)
    setj <- getWithColData(object, name)
    if (any(!subset %in% rowData(seti)[[fcol]]))
        stop("Some entries in 'subset' not found in '", fcol,
             "' (rowData of set '", i, "')")
    setjRefit <- msqrobAggregate(
        seti[rowData(seti)[[fcol]] %in% subset, ],
        formula = formula, fcol = fcol, modelColumnName = modelColumnName,
        ...
    )
    rowData(setj)[[modelColumnName]][subset] <-
        rowData(setjRefit)[[modelColumnName]][subset]
    modelsNew <- rowData(setj)[[modelColumnName]]
    hlp <- limma::squeezeVar(
        var = vapply(modelsNew, getVar, numeric(1)),
        df = vapply(modelsNew, getDF, numeric(1))
    )
    for (ii in seq_along(modelsNew)) {
        modelsNew[[ii]]@varPosterior <- as.numeric(hlp$var.post[ii])
        modelsNew[[ii]]@dfPosterior <- as.numeric(hlp$df.prior + getDF(modelsNew[[ii]]))
    }
    rowData(object[[name]])[[modelColumnName]] <- modelsNew
    object
}


msqrob2Workflow <- function(data, model, robust, ridge, formula, i,
                            psmLevel, technicalReplicates) {
    ## Model data
    data <- runMsqrob2(
        data, model, robust, ridge, formula, i, psmLevel,
        technicalReplicates
    )
    if (psmLevel) i <- "proteins_msqrob2"
    ## Hypothesis testing
    L <- generateContrast(ridge)
    data <- hypothesisTest(
        data, i = i, contrast = L, resultsColumnNamePrefix = "test_"
    )
    ## Format results
    out <- inferenceResults(data, i)
    out$Model <- model
    out$technicalReplicates <- technicalReplicates
    gc() ## free up some RAM space
    out
}

runMsqrob2 <- function(data, model, robust, ridge, formula, i,
                       psmLevel, technicalReplicates) {
    ## Get data: with or without technical replication
    if (!technicalReplicates) {
        data <- subsetByColData(data, data$TechRepMixture == 1)
        data <- dropEmptyAssays(data)
    }
    formula <- as.formula(formula)
    ## Perform psm-level or protein-level modelling
    if (psmLevel) {
        data <- msqrobAggregate(
            data, i = i, fcol = "Protein.Accessions",
            formula = formula, ridge = ridge, robust = robust,
            name = "proteins_msqrob2", aggregateFun = colSums
        )
    } else {
        data <- msqrob(
            data, i = i, formula = formula, ridge = ridge,
            robust = robust
        )
    }
    data
}

generateContrast <- function(ridge) {
    paramNames <- c("Condition1","Condition0.667","Condition0.5","Condition0.125")
    combinations <- sapply(
        combn(paramNames, 2, simplify = FALSE),
        function(x) paste(paste(x, collapse = " - "), " = 0")
    )
    L <- makeContrast(combinations, paramNames)
    if (ridge)
        dimnames(L) <- list(gsub("Cond", "ridgeCond", rownames(L)),
                            gsub("Cond", "ridgeCond", colnames(L)))
    L
}

inferenceResults <- function(data, i) {
    rd <- rowData(data[[i]])
    ## Inference results were stored using the "test_" prefix
    tests <- rd[, grep("^test_", colnames(rd))]
    ## A little cleaning
    tests <- lapply(names(tests), function(ii) {
        test <- tests[[ii]]
        test$Comparison <- sub(".*(ridge)?Condition(.*) - (ridge)?Condition(.*)", "\\2 - \\4", ii)
        test$Protein <- rownames(test)
        test
    })
    ## Return a single table for all comparisons
    do.call(rbind, tests)
}

####---- Plotting utils ----####


colours <- c(
    "DEqMS"="#66c2a5",
    "MSstatsTMT" = "#ffb55a",
    "msqrob2_lmm" = "grey80",
    "msqrob2_rlmm" = "#669e55",
    "msqrob2_rilmm" = "#d6c349",
    "msqrob2_rrilmm" ="#2b8cbe" ,
    "msqrob2_rrilmm_mixture" = "#093c57",
    "msqrob2_psm_rrilmm_refit" = "#bd7ebe",
    "msqrob2_psm_lmm" = "grey80",
    "msqrob2_psm_rlmm" = "#669e55",
    "msqrob2_psm_rilmm" = "#d6c349",
    "msqrob2_psm_rrilmm" = "#fd7f6f",
    "msTrawler" = "grey35",
    "msTrawler_fixed" = "grey70",
    "MSstatsTMT_imp_refnorm" = "#ffb55a",
    "MSstatsTMT_imp_norefnorm" = "#de9337",
    "MSstatsTMT_noimp_refnorm" = "#b36d19",
    "MSstatsTMT_noimp_norefnorm" = "#633701",
    "msqrob2tmt_lmm_imp_refnorm" = "#7bccc4",
    "msqrob2tmt_lmm_imp_norefnorm" = "#3e9c92",
    "msqrob2tmt_lmm_noimp_refnorm" = "#10756a",
    "msqrob2tmt_lmm_noimp_norefnorm" = "#014a40"
)

computeFDP <- function(pval, tp) {
    ord <- order(pval)
    fdp <- cumsum(!tp[ord]) / 1:length(tp)
    fdp[order(ord)]
}

computeTPR <- function(pval, tp, nTP = NULL) {
    if (is.null(nTP)) nTP <- sum(tp)
    ord <- order(pval)
    tpr <- cumsum(tp[ord]) / nTP
    tpr[order(ord)]
}

plotTprFdp <- function(df, colours, groupByComparison = FALSE, nTP = NULL) {
    if (groupByComparison) {
        grouping <- function(df) group_by(df, Model, Comparison)
    } else {
        grouping <- function(df) group_by(df, Model)
    }
    df <- grouping(df) |>
        mutate(tpr = computeTPR(pval, Differential, nTP),
               fdp = computeFDP(pval, Differential)) |>
        arrange(fdp)
    ggplot(df) +
        aes(y = fdp,
            x = tpr,
            fill = Model,
            colour = Model) +
        geom_borderline(
            bordercolour = "black",
            borderwidth  = 0.1,
            linejoin = "mitre"
        ) +
        geom_hline(yintercept = 0.05, linetype = 2, color = "grey") +
        geom_point(
            data = grouping(df) |>
                filter(adjPval < 0.05) |>
                slice_max(adjPval) |>
                filter(!duplicated(Model)),
            size = 3, shape = 21, stroke = 1, colour = "black"
        ) +
        labs(y = "False discovery proportion", x = "True positive rate") +
        coord_flip() +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        theme_bw()
}

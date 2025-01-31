library("dplyr")
library("ggplot2")
library("patchwork")
library("UpSetR")
library("gprofiler2")
library("ggborderline")
library("tidyr")

dataDir <- "data/"
plotDir <- "figs/"

source("scripts/utils.R")

####---- Load data ----####

## List of all comparisons considered for this dataset
comparisons <- c(
    "LF - HF (early)",
    "LF - HF (late)",
    "LF - HF (averaged)",
    "Interaction"
)

inference <- list.files(dataDir, "mouse_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    mutate(
        Comparison = case_match(
            Comparison,
            "Short_LF - Short_HF" ~ "LF - HF (early)",
            "Long_LF - Long_HF" ~ "LF - HF (late)",
            "(Short_LF + Long_LF)/2 - (Short_HF + Long_HF)/2" ~ "LF - HF (averaged)",
            "(Short_LF - Short_HF) - (Long_LF - Long_HF)" ~ "Interaction"),
        Comparison = factor(Comparison, levels = comparisons)
    )

####---- Supplementary: P-value distribution under the null ----####

models <- c(
    "msqrob2_rrilmm_mixture",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)

dfPlotS1 <- mutate(inference,
                   Model = factor(Model, levels = models)) |>
    filter(ModelVariable == "Mock",
           Model %in% models)

(figS1 <- ggplot(dfPlotS1) +
        aes(x = pval, color = Model) +
        geom_borderstep(
            stat = "bin", bins =40, position = "identity",
            direction = "vh", bordercolour = "black",
            borderwidth = 0.01
        ) +
        scale_colour_manual(values = colours) +
        coord_cartesian(y = c(0, 200)) +
        xlab("p-value") +
        guides(col = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.position = "bottom"))
ggsave(
    paste0(plotDir, "mouse_supp_null_distribution.png"), plot = figS1,
    width = 6, height = 4
)

####---- Supplementary: Upset plots ----####

models <- c(
    "MSstatsTMT",
    "msqrob2_rrilmm_mixture",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)

significantHits <- filter(inference,
                          ModelVariable == "Condition",
                          !is.na(adjPval),
                          adjPval < 0.05,
                          Model %in% models)
significantProteinList <- lapply(
    split(significantHits, significantHits$Comparison),
    function(x) {
        x <- filter(x, Model %in% models)
        prots <- split(x$Protein, x$Model)
    })
upsetPlots <- lapply(names(significantProteinList), function(i) {
    l <- significantProteinList[[i]]
    upset(fromList(l), order.by = "freq")
})
names(upsetPlots) <- unique(significantHits$Comparison)

png(paste0(plotDir, "mouse_upset_early.png"), res = 400, width = 2200, height = 2200)
upsetPlots$"LF - HF (early)"
dev.off()

png(paste0(plotDir, "mouse_upset_late.png"), res = 400, width = 2200, height = 2200)
upsetPlots$"LF - HF (late)"
dev.off()

png(paste0(plotDir, "mouse_upset_average.png"), res = 400, width = 2200, height = 2200)
upsetPlots$"LF - HF (averaged)"
dev.off()

png(paste0(plotDir, "mouse_upset_interaction.png"), res = 400, width = 2200, height = 2200)
upsetPlots$Interaction
dev.off()

####---- Supplementary: Over-representation analysis ----####

models <- c(
    "MSstatsTMT",
    "msqrob2_rrilmm_mixture",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)
comparison <- "LF - HF (averaged)"

significantHitsi <- filter(significantHits, Comparison == comparison) ## Focus on the average comparison
significantProteinListi <- split(significantHitsi$Protein, significantHitsi$Model)
protsInData <- unique(readRDS(paste0(dataDir, "mouse_model_msqrob2tmt.rds"))$Protein)

goResAll <- lapply(c("mouse", "identified"), function(background) {
    if (background == "identified") {
        custom_bg <- protsInData
    } else {
        custom_bg <- NULL
    }
    goRes <- lapply(significantProteinListi, function(x) {
        goRes <- gost(
            query = x, organism = "mmusculus", sources = "KEGG",
            custom_bg = custom_bg
        )$result
    })
    goRes <- do.call(rbind, goRes)
    goRes$Background <- background
    goRes$Model <- sub("[.]\\d+", "", rownames(goRes))
    goRes$Model <- factor(goRes$Model, levels = models)
    goRes
})

(figS3 <- do.call(rbind, goResAll) |>
        ggplot() +
        aes(x = Model, y = term_name) +
        geom_tile(fill = "grey50", colour = "grey30", linewidth = 0.5) +
        facet_grid(~ Background, labeller = label_both) +
        labs(title = "Over-representation analysis",
             x = "", y = "") +
        scale_x_discrete(limits = models) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
ggsave(
    paste0(plotDir, "mouse_supp_ora.png"), plot = figS3,
    width = 6, height = 9
)

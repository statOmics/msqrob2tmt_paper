
library("dplyr")
library("ggplot2")
library("ggborderline")
library("tidyr")
library("patchwork")

dataDir <- "data/"
plotDir <- "figs/"

source("scripts/utils.R")

####---- Load data ----####

inferenceAll <- list.files(dataDir, "spikein1_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    filter(!is.na(adjPval)) |>
    select(Protein, logFC, pval, adjPval, Comparison, Model, technicalReplicates) |>
    mutate(Differential = grepl("ups", Protein))

## Benchmark on the data set without technical replicates
inference <- filter(inferenceAll, !technicalReplicates) |>
    select(-technicalReplicates)
## Benchmark on the full data set, with technical replicates
inferenceWithRep <- filter(inferenceAll, technicalReplicates) |>
    select(-technicalReplicates)

####---- Figure 2: TPR-FDP curves ----####

models <- c(
    "DEqMS",
    "MSstatsTMT",
    "msTrawler",
    "msTrawler_fixed",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)

dfPlot2 <- filter(inference, Model %in% models) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           Model = factor(Model, levels = models)) |>
    ungroup()

(fig2A <- plotTprFdp(dfPlot2, colours) +
        ylim(c(0, 0.25)) +
        guides(fill = "none",
               colour = "none"))
(fig2B <- group_by(dfPlot2, Model) |>
        summarise(n = length(unique(Protein))) |>
        ggplot() +
        aes(y = n,
            x = Model,
            fill = Model,
            label = n) +
        geom_bar(stat = "identity") +
        geom_text(vjust = 1.5, size = 2.5, colour = "white") +
        labs(title = "Fitted proteins", x = "", y = "Count") +
        scale_fill_manual(values = colours) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()))

dfPlot2C <- group_by(dfPlot2, Model) |>
    summarise("True Positives" = sum(adjPval < 0.05 & Differential),
              "False Positives" = sum(adjPval < 0.05 & !Differential)) |>
    pivot_longer(cols = c("True Positives", "False Positives"))
(fig2C <- ggplot(dfPlot2C) +
        aes(x = Model,
            y = value,
            fill = Model,
            label = value) +
        geom_bar(stat = "identity") +
        geom_text(data = filter(dfPlot2C, value > 4),
                  vjust = 1.5, size = 2.5, colour = "white") +
        geom_text(data = filter(dfPlot2C, value < 4),
                  vjust = -0.5, size = 2.5) +
        facet_grid(name ~ ., scales = "free") +
        scale_fill_manual(values = colours) +
        labs(title = "5% FDR threshold", x = "", y = "Count") +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()))
(fig2 <- fig2A +
        fig2B +
        fig2C +
        plot_layout(design = "AAB
                              AAC
                              AAC",
                    guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "figure2.png"), plot = fig2,
    width = 8, height = 6
)

####---- Figure 3: LogFC boxplots ----####

dfPlot3 <- mutate(inference,
                  Model = factor(Model, levels = models)) |>
    filter(Model %in% models,
           Comparison %in% c("1 - 0.125", "0.667 - 0.5"))

## Some outliers were found, here is an overview of excluded entries
limMin <- -3
(outlierDf <- filter(dfPlot3, logFC < limMin))

(fig3A <- filter(dfPlot3, !Differential) |>
        ggplot() +
        aes(y = logFC,
            colour = Model) +
        geom_boxplot() +
        geom_hline(yintercept = 0, colour = "black", linetype = 2,
                   linewidth = 0.8, alpha = 0.5) +
        facet_grid(~ Comparison, scales= "free") +
        ylim(limMin, NA) +
        labs(title = "Yeast background proteins", y = "log2 fold change") +
        scale_colour_manual(values = colours) +
        theme_bw())

(fig3B <- filter(dfPlot3, Differential) |>
        ggplot() +
        aes(y = logFC,
            colour = Model) +
        geom_boxplot() +
        geom_hline(
            data = data.frame(
                Comparison = c("0.667 - 0.5","1 - 0.125"),
                logFC = c(log2(0.667 / 0.5), log2(1 / 0.125))
            ),
            aes(yintercept = logFC),
            colour = "black", linetype = 2, linewidth = 0.8,
            alpha = 0.5
        ) +
        labs(title = "UPS spike-in proteins", y = "log2 fold change") +
        facet_wrap(~ Comparison) +
        scale_colour_manual(values = colours) +
        theme_bw())

(fig3 <- fig3A +
        fig3B +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") &
        guides(col = guide_legend(nrow = 2, byrow = TRUE)) &
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank())
)
ggsave(
    paste0(plotDir, "figure3.png"), plot = fig3,
    width = 7, height = 5
)

####---- Figure 4: Compare msqrob2TMT flavours ----####

models <- c(
    "msqrob2_psm_lmm",
    "msqrob2_psm_rlmm",
    "msqrob2_psm_rilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_lmm",
    "msqrob2_rlmm",
    "msqrob2_rilmm",
    "msqrob2_rrilmm"
)
dfPlot4 <- filter(inference, Model %in% models) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           level = ifelse(grepl("psm", Model), "PSM-level models", "Protein-level models"),
           Model = factor(Model, levels = models),
           level = factor(level, levels = c("PSM-level models", "Protein-level models"))) |>
    ungroup()
fig4 <- lapply(split(dfPlot4, dfPlot4$level), function(x) {
    plotTprFdp(x, colours) +
        xlim(c(0.6, 1)) +
        ylim(c(0, 0.25)) +
        facet_grid(~ level) +
        theme(legend.position = "bottom",
              legend.title = element_blank()) +
        guides(colour = guide_legend(nrow = 2))
})
(fig4 <- wrap_plots(fig4))
ggsave(
    paste0(plotDir, "figure4.png"), plot = fig4,
    width = 7, height = 4
)

####---- Figure S1: TPR-FDP curves for 3 strategies ----####

models <- c(
    "DEqMS",
    "MSstatsTMT",
    "msTrawler",
    "msTrawler_fixed",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)

## figS1A: Include only the proteins fit by each model, irrespective
## of the set of proteins fit by the other models
dfPlotS1 <- dfPlot2
(figS1A <- plotTprFdp(dfPlotS1, colours))

## figS1B: Include all ground-truth DA proteins
## Retrieve all the proteins, whether they are fit or not
allProts <- list.files(dataDir, "spikein1_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    pull(Protein) |>
    unique()
## Use the number of nDA as the known number of true positives
nDA <- sum(grepl("ups", allProts)) * length(unique(dfPlotS1$Comparison))
(figS1B <- plotTprFdp(dfPlotS1, colours, nTP = nDA))


## figS1C: Include only the set of proteins found across all models
commonProteins <- split(dfPlotS1$Protein, paste0(dfPlotS1$Model)) |>
    Reduce(f = intersect)
(figS1C <- filter(dfPlotS1, Protein %in% commonProteins) |>
        plotTprFdp(colours))

## Combine panels
(figS1 <- figS1A +
        ggtitle("Performance considering only proteins fit by each method individually") +
        figS1B +
        ggtitle("Performance considering all ground truth DA proteins") +
        figS1C +
        ggtitle("Performance considering proteins commonly fit by all methods") +
        plot_layout(guides = "collect", ncol = 1) +
        plot_annotation(tag_levels = "A") &
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)) &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "figureS1.png"), plot = figS1,
    width = 8, height = 10
)


####---- Figure S2-4: TPR-FDP curves - individual comparisons ----####

dfPlotS2 <- filter(dfPlotS1, Model %in% models) |>
    group_by(Model, Comparison) |>
    mutate(adjPval = p.adjust(pval, method = "BH"))
(figS2 <- plotTprFdp(dfPlotS2, colours, groupByComparison = TRUE) +
        facet_wrap(~ Comparison) +
        theme(legend.position = "bottom",
              legend.title = element_blank()) +
    guides(colour =guide_legend(nrow = 2, byrow = TRUE)))
ggsave(
    paste0(plotDir, "figureS2.png"), plot = figS2,
    width = 8, height = 8
)

nDA <- sum(grepl("ups", allProts))
dfPlotS3 <- filter(dfPlotS1, Model %in% models) |>
    group_by(Model, Comparison) |>
    mutate(adjPval = p.adjust(pval, method = "BH"))
(figS3 <- plotTprFdp(dfPlotS3, colours, groupByComparison = TRUE, nTP = nDA) +
        facet_wrap(~ Comparison) +
        theme(legend.position = "bottom",
              legend.title = element_blank()) +
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)))
ggsave(
    paste0(plotDir, "figureS3.png"), plot = figS3,
    width = 8, height = 8
)

dfPlotS4 <- filter(dfPlotS1,
                   Protein %in% commonProteins,
                   Model %in% models) |>
    group_by(Model, Comparison) |>
    mutate(adjPval = p.adjust(pval, method = "BH"))
(figS4 <- plotTprFdp(dfPlotS4, colours, groupByComparison = TRUE) +
        facet_wrap(~ Comparison) +
        theme(legend.position = "bottom",
              legend.title = element_blank()) +
        guides(colour =guide_legend(nrow = 2, byrow = TRUE)))
ggsave(
    paste0(plotDir, "figureS4.png"), plot = figS4,
    width = 8, height = 8
)

####---- Figure S5: spike-in with replication  ----####

models <- c(
    "DEqMS",
    "MSstatsTMT",
    "msTrawler",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)

## figS5A: Include only the proteins fit by each model, irrespective
## of the set of proteins fit by the other models
dfPlotS5 <- filter(inferenceWithRep, Model %in% models) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           Model = factor(Model, levels = models)) |>
    ungroup()

(figS5A <- plotTprFdp(dfPlotS5, colours))

## figS5B: Include all ground-truth DA proteins
## Retrieve all the proteins, whether they are fit or not
allProts <- list.files(dataDir, "spikein1_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    pull(Protein) |>
    unique()
## Use the number of nDA as the known number of true positives
nDA <- sum(grepl("ups", allProts)) * length(unique(dfPlotS5$Comparison))
(figS5B <- plotTprFdp(dfPlotS5, colours, nTP = nDA))

## figS5C: Include only the set of proteins found across all models
commonProteins <- split(dfPlotS5$Protein, paste0(dfPlotS5$Model)) |>
    Reduce(f = intersect)
(figS5C <- filter(dfPlotS5, Protein %in% commonProteins) |>
        plotTprFdp(colours))

## Combine panels
(figS5 <- figS5A +
        ggtitle("Performance considering only proteins fit by each method individually") +
        figS5B +
        ggtitle("Performance considering all proteins, p-values for non fit proteins set to 1") +
        figS5C +
        ggtitle("Performance considering proteins commonly fit by all methods") +
        plot_layout(guides = "collect", ncol = 1) +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "figureS5.png"), plot = figS5,
    width = 8, height = 10
)

####---- Figure S6: P-value distribution under the null ----####

dfPlotS6 <- mutate(inference,
                  Model = factor(Model, levels = models)) |>
    filter(!Differential,
           Model %in% models,
           Comparison %in% c("1 - 0.125", "0.667 - 0.5"))

(figS6 <- ggplot(dfPlotS6) +
        aes(x = pval, color = Model) +
        geom_borderstep(
            stat = "bin", bins =40, position = "identity",
            direction = "vh", bordercolour = "black",
            borderwidth = 0.01
        ) +
        scale_colour_manual(values = colours) +
        facet_grid(Comparison ~ .) +
        coord_cartesian(y = c(0, 200)) +
        xlab("p-value") +
        guides(col = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.position = "bottom"))
ggsave(
    paste0(plotDir, "figureS6.png"), plot = figS6,
    width = 7, height = 5
)

####---- Figure S7: compare preprocessing workflows ----####

models <- c(
    "msqrob2tmt_lmm_noimp_norefnorm",
    "msqrob2tmt_lmm_noimp_refnorm",
    "msqrob2tmt_lmm_imp_norefnorm",
    "msqrob2tmt_lmm_imp_refnorm",
    "msqrob2_rrilmm",
    "MSstatsTMT_noimp_norefnorm",
    "MSstatsTMT_noimp_refnorm",
    "MSstatsTMT_imp_norefnorm",
    "MSstatsTMT_imp_refnorm"
)

newLabels <- c(
    msqrob2tmt_lmm_imp_refnorm = "lmm_imp_refnorm",
    msqrob2tmt_lmm_imp_norefnorm = "lmm_imp_norefnorm",
    msqrob2tmt_lmm_noimp_refnorm = "lmm_noimp_refnorm",
    msqrob2tmt_lmm_noimp_norefnorm = "lmm_noimp_norefnorm",
    msqrob2_rrilmm = "msqrob2_rrilmm",
    MSstatsTMT_imp_refnorm = "MSstatsTMT_imp_refnorm",
    MSstatsTMT_imp_norefnorm = "MSstatsTMT_imp_norefnorm",
    MSstatsTMT_noimp_refnorm = "MSstatsTMT_noimp_refnorm",
    MSstatsTMT_noimp_norefnorm = "MSstatsTMT_noimp_norefnorm"
)

dfPlotS7 <- filter(inference, Model %in% models) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           Model = factor(Model, levels = models)) |>
    ungroup()

## figS7A: Include only the proteins fit by each model, irrespective
## of the set of proteins fit by the other models
(figS7A <- plotTprFdp(dfPlotS7, colours))

## figS7B: Include all ground-truth DA proteins
## Retrieve all the proteins, whether they are fit or not
allProts <- list.files(dataDir, "spikein1_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    pull(Protein) |>
    unique()
## Use the number of nDA as the known number of true positives
nDA <- sum(grepl("ups", allProts)) * length(unique(dfPlotS7$Comparison))
(figS7B <- plotTprFdp(dfPlotS7, colours, nTP = nDA))

## figS7C: Include only the set of proteins found across all models
commonProteins <- split(dfPlotS7$Protein, paste0(dfPlotS7$Model)) |>
    Reduce(f = intersect)
(figS7C <- filter(dfPlotS7, Protein %in% commonProteins) |>
        plotTprFdp(colours))

## Combine panels
(figS7 <- figS7A +
        ggtitle("Performance considering only proteins\nfit by each method individually") +
        figS7B +
        ggtitle("Performance considering all ground truth DA proteins") +
        figS7C +
        ggtitle("Performance considering proteins\ncommonly fit by all methods") +
        plot_layout(guides = "collect", ncol = 1) +
        plot_annotation(tag_levels = "A") &
        scale_colour_manual(labels = newLabels, values = colours) &
        scale_fill_manual(labels = newLabels, values = colours) &
        ylim(0, 0.25) &
        xlim(0.5, 1) &
        guides(colour = guide_legend(ncol = 2, byrow = FALSE)) &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "figureS7.png"), plot = figS7,
    width = 5.5, height = 10
)


library("dplyr")
library("ggplot2")
library("ggborderline")
library("tidyr")
library("patchwork")

dataDir <- "data/"
plotDir <- "figs/"

source("scripts/utils.R")

####---- Load data ----####

inference <- list.files(dataDir, "spikein2_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    filter(!is.na(adjPval),
           ## Keep only models that were run on their matched preprocessing workflow
           !(grepl("msqrob", Model) & Preprocessing == "msTrawler"))  |>
    select(Protein, logFC, pval, adjPval, Comparison, Model) |>
    mutate(Differential = grepl("YEAST", Protein))

####---- TPR-FDP curves ----####

models <- c(
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit",
    "DEqMS",
    "msTrawler"
)

dfPlot1 <- filter(inference, Model %in% models) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           Model = factor(Model, levels = models)) |>
    ungroup()

(fig1A <- plotTprFdp(dfPlot1, colours) +
        ylim(c(0, 0.25)) +
        guides(fill = "none",
               colour = "none"))
(fig1B <- group_by(dfPlot1, Model) |>
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
dfPlot1C <- group_by(dfPlot1, Model) |>
    summarise("True Positives" = sum(adjPval < 0.05 & Differential),
              "False Positives" = sum(adjPval < 0.05 & !Differential)) |>
    pivot_longer(cols = c("True Positives", "False Positives"))
(fig1C <- ggplot(dfPlot1C) +
        aes(x = Model,
            y = value,
            fill = Model,
            label = value) +
        geom_bar(stat = "identity") +
        geom_text(data = filter(dfPlot1C, value > 200),
                  vjust = 1.5, size = 2.5, colour = "white") +
        geom_text(data = filter(dfPlot1C, value < 200),
                  vjust = -0.5, size = 2.5) +
        facet_grid(name ~ ., scales = "free") +
        scale_fill_manual(values = colours) +
        labs(title = "5% FDR threshold", x = "", y = "Count") +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()))
(fig1 <- fig1A +
        fig1B +
        fig1C +
        plot_layout(design = "AAB
                              AAC
                              AAC",
                    guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "spikein2_tpr_fdp.png"), plot = fig1,
    width = 8, height = 6
)

####---- LogFC boxplots ----####

targetComparisons <- c("11 - 1", "2 - 1")
dfPlot2 <- mutate(inference,
                  Model = factor(Model, levels = models)) |>
    filter(Model %in% models,
           Comparison %in% targetComparisons)

## Some outliers were found, here is an overview of excluded entries
limMax <- 2
(outlierDf <- filter(dfPlot2, logFC > limMax & !Differential))

(fig2A <- filter(dfPlot2, !Differential) |>
        ggplot() +
        aes(y = logFC,
            colour = Model) +
        geom_boxplot() +
        geom_hline(yintercept = 0, colour = "black", linetype = 2,
                   linewidth = 0.8, alpha = 0.5) +
        labs(title = "Mouse background proteins", y = "log2 fold change") +
        facet_grid(~ Comparison, scales= "free") +
        ylim(NA, limMax) +
        scale_colour_manual(values = colours) +
        theme_bw())

(fig2B <- filter(dfPlot2, Differential) |>
        ggplot() +
        aes(y = logFC,
            colour = Model) +
        geom_boxplot() +
        geom_hline(
            data = data.frame(
                Comparison = targetComparisons,
                logFC = c(log2(11), log2(2))
            ),
            aes(yintercept = logFC),
            colour = "black", linetype = 2, linewidth = 0.8,
            alpha = 0.5
        ) +
        labs(title = "Yeast spike-in proteins", y = "log2 fold change") +
        facet_wrap(~ Comparison) +
        scale_colour_manual(values = colours) +
        theme_bw())

(fig2 <- fig2A +
        fig2B +
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
    paste0(plotDir, "spikein2_logFC.png"), plot = fig2,
    width = 7, height = 5
)

####---- P-value distribution under the null ----####

dfPlot3 <- mutate(inference,
                  Model = factor(Model, levels = models)) |>
    filter(!Differential,
           Model %in% models,
           Comparison %in% c("11 - 1", "2 - 1", "14 - 11", "20 - 18"))

(fig3 <- ggplot(dfPlot3) +
        aes(x = pval, color = Model) +
        geom_borderstep(
            stat = "bin", bins =40, position = "identity",
            direction = "vh", bordercolour = "black",
            borderwidth = 0.01
        ) +
        scale_colour_manual(values = colours) +
        facet_grid(Comparison ~ .) +
        coord_cartesian(y = c(0, 20)) +
        xlab("p-value") +
        guides(col = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.position = "bottom"))
ggsave(
    paste0(plotDir, "spikein2_null_distribution.png"), plot = fig3,
    width = 7, height = 5
)

####---- Supplementary: TPR-FDP curves for 3 strategies ----####

models <- c(
    "DEqMS",
    "MSstatsTMT",
    "msTrawler",
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit"
)

## figS1A: Include only the proteins fit by each model, irrespective
## of the set of proteins fit by the other models
dfPlotS1A <- filter(inference, Model %in% models) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           Model = factor(Model, levels = models)) |>
    ungroup()
(figS1A <- plotTprFdp(dfPlotS1A, colours))

## figS1B: Include all ground-truth DA proteins
## Retrieve all the proteins, whether they are fit or not
allProts <- list.files(dataDir, "spikein2_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    pull(Protein) |>
    unique()
## Use the number of nDA as the known number of true positives
nDA <- sum(grepl("YEAST", allProts)) * length(unique(dfPlotS1A$Comparison))
(figS1B <- plotTprFdp(dfPlotS1A, colours, nTP = nDA))

## figS1C: Include only the set of proteins found across all models
commonProteins <- split(dfPlotS1A$Protein, paste0(dfPlotS1A$Model)) |>
    Reduce(f = intersect)
(figS1C <- filter(dfPlotS1A, Protein %in% commonProteins) |>
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
        ylim(0, 0.25) &
        xlim(0, 1) &
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)) &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "spikein2_supp_tpr_fdp.png"), plot = figS1,
    width = 8, height = 10
)

####---- Supplementary: compare preprocessing workflows ----####

models <- c(
    "msqrob2_rrilmm",
    "msqrob2_psm_rrilmm",
    "msqrob2_psm_rrilmm_refit",
    "msTrawler"
)

dfPlotS2 <- list.files(dataDir, "spikein2_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    filter(!is.na(adjPval),
           Preprocessing == "msTrawler",
           Model %in% models) |>
    select(Protein, logFC, pval, adjPval, Comparison, Model) |>
    mutate(Differential = grepl("YEAST", Protein)) |>
    group_by(Model) |>
    mutate(adjPval = p.adjust(pval, method = "BH"),
           Model = factor(Model, levels = models)) |>
    ungroup()

## figS2A: Include only the proteins fit by each model, irrespective
## of the set of proteins fit by the other models
(figS2A <- plotTprFdp(dfPlotS2, colours))

## figS2B: Include all ground-truth DA proteins
## Retrieve all the proteins, whether they are fit or not
allProts <- list.files(dataDir, "spikein2_model", full.names = TRUE) |>
    lapply(readRDS) |>
    do.call(what = rbind) |>
    pull(Protein) |>
    unique()
## Use the number of nDA as the known number of true positives
nDA <- sum(grepl("YEAST", allProts)) * length(unique(dfPlotS2$Comparison))
(figS2B <- plotTprFdp(dfPlotS2, colours, nTP = nDA))

## figS1C: Include only the set of proteins found across all models
commonProteins <- split(dfPlotS2$Protein, paste0(dfPlotS2$Model)) |>
    Reduce(f = intersect)
(figS2C <- filter(dfPlotS2, Protein %in% commonProteins) |>
        plotTprFdp(colours))

## Combine panels
(figS2 <- figS2A +
        ggtitle("Performance considering only proteins fit by each method individually") +
        figS2B +
        ggtitle("Performance considering all ground truth DA proteins") +
        figS2C +
        ggtitle("Performance considering proteins commonly fit by all methods") +
        plot_layout(guides = "collect", ncol = 1) +
        plot_annotation(tag_levels = "A",
                        title = "Comparing performance after msTrawler preprocessing",
                        theme = theme(title = element_text(size = 12, face = "bold"))) &
        ylim(0, 0.25) &
        xlim(0, 1) &
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)) &
        theme(legend.position = "bottom",
              legend.title = element_blank()))
ggsave(
    paste0(plotDir, "spikein2_supp_compare_preprocess.png"), plot = figS2,
    width = 8, height = 10
)


## Pre-process the spikein2 data set ##

## This script loads the yeast spike-in data by O'Brien et al. 2020
## and performs a custom data pre-processing workflow that accounts
## for the peculiarities of this data set. The code to perform
## preprocessing is provided from the article's repo:
## https://console.cloud.google.com/storage/browser/mstrawler_paper/Code_for_Drive/Interplex_Analyses?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&inv=1&invt=AboFkQ
## These data are saved and will be used by all methods for comparing
## performance.
##
## We will then proceed with the preprocessing performed within
## `msTrawl()`. This further preprocessed data will be analysed by
## msqrob2, allowing to decouple changes in performance associated
## with preprocessing and with modelling between msqrob2 and
## msTrawler.

####---- Setup environment ----####

library("tidyverse")
library("msTrawler")
library("BiocFileCache")

dataDir <- "data/"

####---- Load data ----####

DF <- BiocFileCache() |>
    bfcrpath("https://zenodo.org/records/14767905/files/spikein2_psms.csv?download=1") |>
    read.csv(stringsAsFactors = FALSE)

####---- O'Brien's custom data preprocessing ----####

## The remaineder of the script is mostly copy-pasted from the
## original paper's repository:
## https://console.cloud.google.com/storage/browser/mstrawler_paper.

##oooooooooooooooooooooooooooooooooooooo##
## START COPY-PASTE FROM O'BRIEN'S CODE ##
##oooooooooooooooooooooooooooooooooooooo##

#Clean up the data prior to analysis
#Remove reverse hits and contaminants
revI <- grep("##" , DF$Protein.ID)
contamI <- grep("contam", DF$Protein.ID)
cleanPep <- DF[-union(revI, contamI), ]

#Remove peptidesd to prevent confusion about species labels)
#First remove peptides with ambiguous protein assignment
redunI <- which(cleanPep$Redun > 0)
cleanPep <- cleanPep[-redunI, ]

#Now look for mouse peptides that correlate with the yeast diluton profile
iMat <- cleanPep[ , grep("Adjusted", colnames(cleanPep))]
snMat <- cleanPep[ , grep("Sn", colnames(cleanPep))]

ratiosY <- matrix(c(18,32,11,28,11,11,14,28,24,2,32,20,18,1,1,
                    20,3,3,14,2,24,14,2,32,5,28,3,1,1,1,
                    32,1,11,14,1,3,5,28,5,5,32,11,18,1,1,
                    14,11,20,18,18,20,11,0,24,2,24,28,24,1,1,
                    20,24,20,3,1,2,18,2,5,32,3,3,18,1,1,
                    5,20,5,32,14,24,28,14,1,2,28,1,1,1,1), nrow=6, byrow = T)

uPlex <- c("plex-1", "plex-2", "plex-3", "plex-4", "plex-5", "plex-6")
yCorr <- sapply(1:nrow(iMat), function(x) cor(as.numeric(iMat[x, 1:13]),
                                              ratiosY[match(cleanPep$Plex[x], uPlex), 1:13],
                                              use = "na.or.complete"))

highCorMouse <- which(cleanPep$Species == "Mouse" & yCorr > 0.25)
cleanPep <- cleanPep[-highCorMouse, ]
iMat <- iMat[-highCorMouse, ]
snMat <- snMat[-highCorMouse, ]

#Do two sets of signal noise filtering.
#We are making 2 datasets to evaluate our algorithms
ssnI <- which(apply(snMat, 1, sum) < 20)
cleanPep <- cleanPep[-ssnI, ]
iMat <- iMat[-ssnI, ]
snMat <- snMat[-ssnI, ]

ssnI <- which(apply(snMat, 1, sum) < 200)
cleanPep2 <- cleanPep[-ssnI, ]
iMat2 <- iMat[-ssnI, ]
snMat2 <- snMat[-ssnI, ]

#Now do column adjustments
DF <- cleanPep

speciesI <- which(DF$Species == "Mouse")
pepSD <- apply(iMat, 1, function(x)
    sd(log(x[-15]) - log(x[15]), na.rm=T))
medSD <- median(pepSD[speciesI], na.rm = T)
sdI <- which(pepSD < medSD)
normIndex <- intersect(sdI, speciesI)

normBool <- rep(0, nrow(DF))
normBool[normIndex] <- 1

ratios <- c(rep(1, ncol(snMat) - 1), 1)

normed <- msTrawler::geoNorm(iMat, normBool, DF$Plex, ratios)


DF <- cleanPep2
speciesI <- which(DF$Species == "Mouse")
pepSD <- apply(iMat2, 1, function(x)
    sd(log(x[-15]) - log(x[15]), na.rm=T))
medSD <- median(pepSD[speciesI], na.rm = T)
sdI <- which(pepSD < medSD)
normIndex <- intersect(sdI, speciesI)

normBool <- rep(0, nrow(DF))
normBool[normIndex] <- 1
ratios <- c(rep(1, ncol(snMat2) - 1), 1)

normed2 <- msTrawler::geoNorm(iMat2, normBool, DF$Plex, ratios)

DF1 <- cleanPep
DF1[ , grep("Adjusted", colnames(DF1))] <- normed[[1]]
#write.csv(DF1, "AdjustedPeptides_SN20.csv", row.names = FALSE)


#DF200 <- cleanPep2
#DF200[ , grep("Adjusted", colnames(DF200))] <- normed2[[1]]
#write.csv(DF200, "AdjustedPeptides_SN200.csv", row.names = FALSE)

#Write one more file to compare with MSstats which does not allow repeat peptides
DFu <- DF1
DFu$pepID <- paste0(DFu$Peptide, "_", DFu$Plex)

SSI <- apply(DFu[ , grep("Adjusted", colnames(DFu))], 1, sum, na.rm=T)
DFu <- DFu[order(SSI, decreasing = T) , ]
pepVec <- ""
keepVec <- rep(0, nrow(DFu))

for(i in 1:nrow(DFu)){
    if(DFu$pepID[i] %in% pepVec){
        next
    }else{
        pepVec <- c(pepVec, DFu$pepID[i])
        keepVec[i] <- 1
    }
}

DFu <- DFu[which(keepVec == 1), ]
DFu <- DFu[ , -grep("pepID", colnames(DFu))]

DF <- DFu

##oooooooooooooooooooooooooooooooooooo##
## END COPY-PASTE FROM O'BRIEN'S CODE ##
##oooooooooooooooooooooooooooooooooooo##

## This will be used by all methods as input
saveRDS(DF, paste0(dataDir, "spikein2_input.rds"))

####---- msTrawler's data preprocessing ----####

## In `msTrawl()`, further preprocessing is performed. In this
## section, we copy pasted the preprocessing part of the function to
## reproduce their workflow.

## Function arguments
scaleSN = 1 ## default
lod = 0.01 ## default
imputePenalty = 1 ## default
minAbove = 3 ## as in author's script for multibatch (spikein2) analysis
ssnFilter = NULL ## disable, already performed
outlierCutoff = NULL ## disable, already performed
N_SUM = 3 ## default
swapProtein = FALSE ## default
maxPep = 25 ## default
colAdjust = NULL ## disable, already performed
colRatios = NULL ## disable, already performed
dropContam = FALSE ## disable, already performed
dropReverse = FALSE ## disable, already performed
peptideAnalysis = FALSE ## default
minRE = 5 ## default
timeDiff = FALSE ## not applicable

##ooooooooooooooooooooooooooooooooo##
## START COPY-PASTE FROM MSTRAWL() ##
##ooooooooooooooooooooooooooooooooo##

if (!is.null(covariateFile)) {
    covariateFile <- covariateFile[order(match(covariateFile$Covariate,
                                               colnames(sampleFile))), ]
    coVector <- covariateFile$Covariate
    sampVector <- colnames(sampleFile)
    nMatches <- sum(coVector %in% sampVector)
    if (nMatches < length(coVector)) {
        stop("Error: At least one of the covariate names does not match across files.")
    }
}
if (!is.null(covariateFile)) {
    bridgeI <- grep("BRIDGE", toupper(colnames(sampleFile)))
    if (length(bridgeI) > 0) {
        nMissing <- sum(is.na(sampleFile[-which(sampleFile[,
                                                           bridgeI] == 1), ]))
    }
    else {
        nMissing <- sum(is.na(sampleFile))
    }
    if (nMissing > 0) {
        stop("Missing values are not allowed in covariates")
    }
}
if (length(colAdjust) > 1) {
    DF$colAdjust <- colAdjust
}
if (dropReverse == TRUE) {
    revI <- grep("##", DF$Protein.ID)
    if (length(revI) > 0) {
        DF <- DF[-revI, ]
    }
}
if (dropContam == TRUE) {
    contamI <- grep("contam", DF$Protein.ID)
    if (length(contamI) > 0) {
        DF <- DF[-contamI, ]
    }
}
DF <- DF[order(DF$Protein.ID, DF$Plex, DF$Peptide), ]
if (peptideAnalysis) {
    DF$PA.Gene.Symbol <- paste0(DF$PA.Gene.Symbol, "___",
                                DF$Protein.ID)
    DF$Protein.ID <- DF$Peptide
}
snMat <- DF[, grep("\\.Sn", colnames(DF))]
iMat <- DF[, grep("Adjusted", colnames(DF))]
if (ncol(snMat) != ncol(iMat)) {
    stop("Error: The number of SNR columns\n    does not equal the number of intensity columns. Often this occurs\n    because of an unwanted column containing the string \".Sn\" in the\n                                     column name.")
}
if (!is.null(ssnFilter)) {
    lowI <- which(apply(snMat, 1, sum) < ssnFilter)
    if (length(lowI) > 0) {
        DF <- DF[-lowI, ]
        snMat <- snMat[-lowI, ]
        iMat <- iMat[-lowI, ]
    }
}
emptyRow <- which(apply(iMat, 1, function(x) {
    sum(x > 0, na.rm = T)
}) < 2)
if (length(emptyRow > 0)) {
    DF <- DF[-emptyRow, ]
    iMat <- iMat[-emptyRow, ]
    snMat <- snMat[-emptyRow, ]
}
lodRES <- tmtLOD(snMat, iMat, lod, minAbove, scaleSN, imputePenalty)
tooFew <- lodRES[[3]]
iMat <- lodRES[[1]]
snMat <- lodRES[[2]]
if (sum(tooFew) > 0) {
    DF <- DF[-which(tooFew == TRUE), ]
    iMat <- iMat[-which(tooFew == TRUE), ]
    snMat <- snMat[-which(tooFew == TRUE), ]
}
DF$Protein.ID <- as.character(DF$Protein.ID)
DF$Peptide <- as.character(DF$Peptide)
uPlex <- unique(DF$Plex)
uProt <- unique(DF$Protein.ID)
for (j in 1:length(uPlex)) {
    for (i in 1:length(uProt)) {
        protIndex <- which(DF$Protein.ID == uProt[i] & DF$Plex ==
                               uPlex[j])
        if (length(protIndex) == 0) {
            next
        }
        subDat <- iMat[protIndex, , drop = F]
        subSn <- snMat[protIndex, , drop = F]
        subPep <- DF[protIndex, "Peptide"]
        subSSN <- apply(subSn, 1, sum)
        if (!is.null(outlierCutoff)) {
            outlierMat <- findOutliers(subDat, subSn, outlierCutoff,
                                       scaleSN)
        }
        else {
            outlierMat <- matrix(0, nrow = nrow(subDat),
                                 ncol = ncol(subDat))
        }
        outlierRows <- apply(outlierMat, 1, sum, na.rm = T)
        outlierIndex <- protIndex[which(outlierRows > 0)]
        if (length(outlierIndex) > 0) {
            if (swapProtein == TRUE) {
                DF[outlierIndex, "Protein.ID"] <- paste0(DF[outlierIndex,
                                                            "Protein.ID"], "___", DF[outlierIndex, "Peptide"])
            }
            else {
                DF[outlierIndex, "Protein.ID"] <- "OUTLIER_REMOVE_AT_ONCE!"
            }
        }
        n_after <- length(which(outlierRows == 0))
        if ((n_after < N_SUM & n_after > 1) | (peptideAnalysis ==
                                               TRUE & n_after > 1)) {
            newSnVec <- apply(subSn[which(outlierRows ==
                                              0), , drop = FALSE], 2, sum)
            newFluxVec <- apply(subDat[which(outlierRows ==
                                                 0), , drop = FALSE], 2, sum)
            keepIndex <- protIndex[which(outlierRows ==
                                             0)][1]
            DF[keepIndex, "Peptide"] <- paste0(DF[protIndex[1],
                                                  "Peptide"], "_SUM")
            DF[keepIndex, grep("\\.Sn", colnames(DF))] <- newSnVec
            DF[keepIndex, grep("Adjusted", colnames(DF))] <- newFluxVec
            snMat[keepIndex, ] <- newSnVec
            iMat[keepIndex, ] <- newFluxVec
            DF[setdiff(protIndex, keepIndex), "Protein.ID"] <- "OUTLIER_REMOVE_AT_ONCE!"
            next
        }
        if (nrow(subDat) - length(outlierIndex) > maxPep) {
            listOranks <- lapply(unique(subPep), function(x) {
                order(subSSN[which(subPep == x)], decreasing = TRUE)
            })
            rankVector <- do.call(c, listOranks)
            if (length(outlierIndex) > 0) {
                rankVector[which(outlierRows > 0)] <- 999999
            }
            snOrder <- order(rankVector, (-1 * subSSN))
            bigIndex <- head(snOrder, n = maxPep)
            smallIndex <- protIndex[setdiff(snOrder, bigIndex)]
            DF[smallIndex, "Protein.ID"] <- "TOO_MANY!!!"
        }
    }
}
outlierIndex <- grep("OUTLIER_REMOVE_AT_ONCE!", DF$Protein.ID)
tooManyIndex <- grep("TOO_MANY!!!", DF$Protein.ID)
removeIndex <- union(outlierIndex, tooManyIndex)
if (length(removeIndex > 0)) {
    DF <- DF[-removeIndex, ]
    iMat <- iMat[-removeIndex, ]
    snMat <- snMat[-removeIndex, ]
}
if (is.null(colAdjust)) {
    normI <- iMat
} else {
    if (length(colAdjust) == 1) {
        normBool <- rep(0, nrow(iMat))
        for (p_ in 1:length(uPlex)) {
            plexIndex <- which(DF$Plex == uPlex[p_])
            if (length(plexIndex) <= 1) {
                stop("You have a plex with <= 1 observation")
            }
            pepSD <- apply(iMat[plexIndex, ], 1, function(x) sd(log(x),
                                                                na.rm = TRUE))
            medSD <- quantile(pepSD, probs = colAdjust,
                              na.rm = T)
            sdI <- which(pepSD < medSD)
            normBool[plexIndex[sdI]] <- 1
        }
    }
    else {
        normBool <- DF$colAdjust
    }
    if (length(colRatios) < 2) {
        ratios <- rep(1, length(uPlex) * ncol(iMat))
    }
    normed <- geoNorm(iMat, normBool, Plex = DF$Plex, ratios)
    normI <- normed[[1]]
    colnames(normI) <- paste0("norm_", colnames(normI))
    normFacs <- normed[[2]]
    write.csv(normFacs, file = "ColAdjustmentFactors.csv",
              row.names = FALSE)
}
usedPlexes <- intersect(uPlex, substring(sampleFile$SampleID,
                                         1, regexpr("_", sampleFile$SampleID) - 1))
nBridges <- sum(sampleFile$Bridge)
if (nBridges > 0) {
    bridgeMod <- TRUE
    if (nBridges != length(usedPlexes)) {
        stop("Error: The number of plexes must match the number of bridge\n         samples.")
    }
} else {
    bridgeMod <- FALSE
}
tableNames <- "Simple.csv"
bridgeCol <- grep("BRIDGE", toupper(colnames(sampleFile)))
if (length(bridgeCol) > 0) {
    coVars <- colnames(sampleFile)[-c(1, bridgeCol)]
} else {
    coVars <- colnames(sampleFile)[-1]
}
if (length(coVars) == 0) {
    coVars <- "1"
}
covType <- covariateFile$Type
covLevel <- covariateFile$Levels
catIndex <- grep("FACTOR", toupper(covType))
n_cat <- length(catIndex)
tParm <- ""
if (n_cat > 0) {
    n_levels <- rep(0, length(catIndex))
    factorNames <- as.character(covariateFile$Covariate[catIndex])
    if (sum(factorNames %in% colnames(sampleFile)) != length(factorNames)) {
        stop("Error: The covariate names in the sampleFile do not match the\n         column names in the sampleFile")
    }
    levelNames <- list()
    for (i in 1:n_cat) {
        if (length(bridgeCol) > 0) {
            realIndex <- which(sampleFile[, bridgeCol] ==
                                   0)
        }
        else {
            realIndex <- 1:nrow(sampleFile)
        }
        levelNames[[i]] <- as.character(unique(sampleFile[realIndex,
                                                          factorNames[i]]))
        n_levels[i] <- length(levelNames[[i]])
        tableNames <- c(tableNames, paste0("Factor_", factorNames[i],
                                           "_", levelNames[[i]], ".csv"))
    }
} else {
    n_levels <- 0
    factorNames <- ""
}
contIndex <- grep("CONTINUOUS", toupper(covariateFile$Type))
n_cont <- length(contIndex)
if (n_cont > 0) {
    contNames <- as.character(covariateFile$Covariate[contIndex])
    tableNames <- c(tableNames, paste0("Continuous_", contNames,
                                       ".csv"))
    for (i in 1:n_cont) {
        sampleFile[, contNames[i]] <- sampleFile[, contNames[i]] -
            mean(sampleFile[, contNames[i]], na.rm = T)
    }
} else {
    contNames <- ""
}
timeIndex <- grep("TIME", toupper(covariateFile$Type))
if (length(timeIndex) > 0) {
    colnames(sampleFile)[grep(covariateFile$Covariate[timeIndex],
                              colnames(sampleFile))] <- "Time"
    covariateFile$Covariate[timeIndex] <- "Time"
    if (!is.numeric(sampleFile$Time)) {
        stop("Error:  Your time variable is non-numeric")
    }
    coVars <- coVars[-timeIndex]
    timeVars <- NULL
    timeDegree <- covariateFile$TimeDegree[timeIndex]
    if (timeDegree >= 1) {
        timeVars <- c(timeVars, "Time")
    }
    if (timeDegree >= 2) {
        timeVars <- c(timeVars, "Time2")
        sampleFile$Time2 <- sampleFile$Time^2
    }
    if (timeDegree == 3) {
        timeVars <- c(timeVars, "Time3")
        sampleFile$Time3 <- sampleFile$Time^3
    }
    tCatIndex <- which(covariateFile$TimeCategory == 1)
    tCatFactorIndex <- match(covariateFile$Covariate[tCatIndex],
                             factorNames)
    if (length(tCatIndex) > 1) {
        stop("Only one time category is allowed")
    }
    if (length(tCatIndex) == 0) {
        tCatIndex <- 0
        tCatFactorIndex <- 0
    }
    if (tCatIndex > 0) {
        tCatName <- covariateFile$Covariate[tCatIndex]
        timeLevelN <- n_levels[tCatFactorIndex]
        tCatLevels <- levelNames[[tCatFactorIndex]]
        timeTableNames <- paste0("Time_", tCatLevels, ".csv")
        tParm <- "Category"
    }
    else {
        timeLevelN <- 1
        timeTableNames <- "Time.csv"
        if (n_cont > 0) {
            tParm <- "Continuous"
        }
        else {
            tParm <- "Time"
        }
    }
    circadian <- covariateFile$Circadian[timeIndex]
} else {
    timeDegree <- 0
    timeVars <- NULL
    timeLevelN <- 0
    tCatIndex <- 0
    tCatFactorIndex <- 0
    circadian <- 0
    tParm <- ""
}
if (circadian != 0) {
    sampleFile$Sine <- sin((2 * pi/24) * sampleFile$Time)
    sampleFile$Cosine <- cos((2 * pi/24) * sampleFile$Time)
    timeVars <- c(timeVars, "Sine", "Cosine")
}
IDindex <- grep("ID", toupper(covariateFile$Type))
if (length(IDindex) > 0) {
    IDname <- covariateFile$Covariate[IDindex]
    randID <- IDname
    coVars <- coVars[-grep(IDname, coVars)]
} else {
    randID <- "SampleID"
}
fixedStr <- paste0("lIntensity ~ ", paste(c(coVars, timeVars),
                                          collapse = " + "))
if (tCatIndex > 0) {
    fixedStr <- paste0(fixedStr, " + ", paste0(c(rep(paste0(tCatName,
                                                            ":"), length(timeVars))), timeVars, collapse = " + "))
}
fixedForm <- as.formula(fixedStr)
idVars <- DF[, c("PA.Gene.Symbol", "Protein.ID", "Peptide",
                 "Plex")]
meltI <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars,
                                   normI), id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID",
                                                       "Peptide", "Plex"))
meltSn <- reshape2::melt(data.frame(Scan = 1:nrow(DF), idVars,
                                    snMat), id.vars = c("Scan", "PA.Gene.Symbol", "Protein.ID",
                                                        "Peptide", "Plex"), value.name = "SNR", variable.name = "Channel")
meltSn$lIntensity <- log2(meltI$value)
meltSn$SampleID <- paste0(meltSn$Plex, "_", meltSn$Channel)
covNames <- colnames(sampleFile)[-1]
if (length(covNames) > 0) {
    dataTosampFile <- match(meltSn$SampleID, sampleFile$SampleID)
    undescribed <- which(is.na(dataTosampFile))
    if (length(undescribed) > 0) {
        meltSn <- meltSn[-undescribed, ]
    }
    idMatch <- match(meltSn$SampleID, sampleFile$SampleID)
    sampFileToData <- match(sampleFile$SampleID, meltSn$SampleID)
    if (sum(is.na(sampFileToData)) > 0) {
        stop("Error: Samples listed in the sampleFile are not present in the data")
    }
    emptyCov <- as.data.frame(matrix(NA, nrow = nrow(meltSn),
                                     ncol = length(covNames)))
    colnames(emptyCov) <- covNames
    for (i in 1:length(covNames)) {
        emptyCov[, i] <- sampleFile[idMatch, i + 1]
    }
    readyDf <- data.frame(meltSn, emptyCov)
} else {
    readyDf <- meltSn
}
naIndex <- which(is.na(readyDf$lIntensity))
if (length(naIndex) > 0) {
    readyDf <- readyDf[-naIndex, ]
}
if (bridgeMod == FALSE) {
    scanMean <- tapply(readyDf$lIntensity, readyDf$Scan,
                       FUN = mean)
    scanFactors <- scanMean
    readyDf$lIntensity <- readyDf$lIntensity - scanFactors[match(readyDf$Scan,
                                                                 names(scanFactors))]
    bridgeDat <- NULL
} else {
    bridgeIDs <- sampleFile$SampleID[which(sampleFile$Bridge ==
                                               1)]
    readyDf$Bridge <- 0
    readyDf$Bridge[which(readyDf$SampleID %in% bridgeIDs)] <- 1
}
sampleNames <- levels(factor(readyDf$SampleID))
readyDf$Plex <- factor(readyDf$Plex)
readyDf$techVar <- 1/(scaleSN * readyDf$SNR)

##ooooooooooooooooooooooooooooooo##
## END COPY-PASTE FROM MSTRAWL() ##
##ooooooooooooooooooooooooooooooo##

## This will be used only by msqrob2tmt methods
saveRDS(readyDf, paste0(dataDir, "spikein2_input_preprocessed.rds"))

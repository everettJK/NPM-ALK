library(GenomicRanges)
library(STRINGdb)
library(biomaRt)
library(DESeq2)
library(xlsx)
library(gtools)
library(tidyverse)
library(pathview)
library(ggrepel)
library(RColorBrewer)
library(xlsx)
source('./lib.R')
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE) 
write(c(date(), capture.output(sessionInfo())), file = 'sessionInfo.txt')
load('data/data.RData')

# Create a report list object to store data needed for reports.
report <- list()
report$STRINGdb_dataFiles <- 'data/STRINGdb' # Created on the fly if not defined.


# STRINGdb will be used for pathway and GO term enrichment.
string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory = report$STRINGdb_dataFiles)


# Read in a list of gene aliases created from HGNC to expand STRINGdb's annotations.
# This data file is created in data_import.R.
report$string_db.alt.aliases <- readRDS('data/string_db.alt.aliases.rds')


# Create DEseq2 data objects which can be used to retrieve counts and create contrasts.
salmon1.ddsTxi <- processSalmon(salmon1.txi)
salmon2.ddsTxi <- processSalmon(salmon2.txi)


# Inspect RNAseq data and select genes to see relationships between counts and fold changes / p-values.
# Cutting the data at +/- 15 fold change appears to be a reasonable means of excluding outliers not caught by DESeq2.
plotMA(salmon1.ddsTxi, ylim = c(-15, 15))
plotMA(salmon2.ddsTxi, ylim = c(-15, 15))
plotCounts(salmon2.ddsTxi, gene="ENSG00000003147.18", intgroup="repGrp")  # ICA1    padj 1e-06    fold change:  2
plotCounts(salmon2.ddsTxi, gene="ENSG00000063046.17", intgroup="repGrp")  # EIF4B   padj 2e-20    fold change: -2
plotCounts(salmon2.ddsTxi, gene="ENSG00000165029.15", intgroup="repGrp")  # ABCA1   padj 8e-230   fold change:  7
plotCounts(salmon2.ddsTxi, gene="ENSG00000080503.23", intgroup="repGrp")  # SMARCA2 padj 5e-88    fold change: -2
plotCounts(salmon2.ddsTxi, gene="ENSG00000077080.9",  intgroup="repGrp")  # ACTL6B  padj NA       fold change: 0.9
plotCounts(salmon2.ddsTxi, gene="ENSG00000117713.19", intgroup="repGrp")  # ARID1A  padj 0.001    fold change: -1
plotCounts(salmon1.ddsTxi, gene="ENSG00000277067.4",  intgroup="repGrp")  # ARID1A  padj 2.97e-07 fold change: 39




# RNAseq PCA plots
#--------------------------------------------------------------------------------------------------

# First RNAseq/SALMON run PCA.
salmon1.rld     <- rlog(salmon1.ddsTxi)  # (Slow) Transform count data with regularized logarithm which returns log2 data normalized to library size. 
salmon1.rlogMat <- assay(salmon1.rld)    # assay() extracts the matrix of normalized values.


# Exclude outlier data points which were not caught by SALMON.
# Threshold of 15x determined by reviewing fold change vs. normalized count plots.
i <- unname(apply(salmon1.rlogMat, 1, function(x) all(x < 15))) 
salmon1.pca <- prcomp(t(salmon1.rlogMat[i,]), scale = FALSE, center = TRUE)


# Use the 'genotype_timePoint_donor' formatted data points to extract data for plots.
salmon1.pca.plotData          <- data.frame(s = row.names(salmon1.pca$x), x = salmon1.pca$x[,1], y = salmon1.pca$x[,2], z = salmon1.pca$x[,3])
salmon1.pca.plotData$genotype <- do.call(rbind, strsplit(salmon1.pca.plotData$s, '_'))[,1]
salmon1.pca.plotData$day      <- do.call(rbind, strsplit(salmon1.pca.plotData$s, '_'))[,2]
salmon1.pca.plotData$subject  <- do.call(rbind, strsplit(salmon1.pca.plotData$s, '_'))[,3]
salmon1.pca.plotData$subject  <- gsub('don', '', salmon1.pca.plotData$subject)
salmon1.pca.plotData          <- salmon1.pca.plotData[mixedorder(as.character(salmon1.pca.plotData$day)),]
salmon1.pca.plotData$day      <- factor(salmon1.pca.plotData$day, levels = unique(salmon1.pca.plotData$day))
salmon1.pca.plotData$genotype <- factor(salmon1.pca.plotData$genotype, levels = c('none', 'KD', 'WT', 'Y664F'))
salmon1.pca.plotData$time     <- as.integer(str_extract(salmon1.pca.plotData$day, '\\d+'))

report$salmon1.pca.plot1 <- 
  make_square(ggplot(salmon1.pca.plotData, aes(x=x, y=y, fill = time, shape = genotype, label = subject)) +
  theme_bw() +
  geom_point(size = 4, stroke = 0.5, color = 'black') +
  scale_shape_manual(name = 'Genotype', values = 21:25) +
  scale_fill_gradient(name = 'Time (days)', low = "white", high = "black") +
  labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
       y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))

ggsave(report$salmon1.pca.plot1, file = 'figures_and_tables/RNAseq1_with_Y664F.pdf')



# WT only plot.
salmon1.rlogMat.noY664F <- salmon1.rlogMat[, -grep('Y664F', colnames(salmon1.rlogMat))]
salmon1.noY664F.pca <- prcomp(t(salmon1.rlogMat.noY664F), scale = FALSE, center = TRUE)

salmon1.pca.plotData2          <- data.frame(s = row.names(salmon1.noY664F.pca$x), x = salmon1.noY664F.pca$x[,1], y = salmon1.noY664F.pca$x[,2], z = salmon1.noY664F.pca$x[,3])
salmon1.pca.plotData2$genotype <- do.call(rbind, strsplit(salmon1.pca.plotData2$s, '_'))[,1]
salmon1.pca.plotData2$day      <- do.call(rbind, strsplit(salmon1.pca.plotData2$s, '_'))[,2]
salmon1.pca.plotData2$subject  <- do.call(rbind, strsplit(salmon1.pca.plotData2$s, '_'))[,3]
salmon1.pca.plotData2$subject  <- gsub('don', '', salmon1.pca.plotData2$subject)
salmon1.pca.plotData2          <- salmon1.pca.plotData2[mixedorder(as.character(salmon1.pca.plotData2$day)),]
salmon1.pca.plotData2$day      <- factor(salmon1.pca.plotData2$day, levels = unique(salmon1.pca.plotData2$day))
salmon1.pca.plotData2$genotype <- factor(salmon1.pca.plotData2$genotype, levels = c('none', 'KD', 'WT'))
salmon1.pca.plotData2$time     <- as.integer(str_extract(salmon1.pca.plotData2$day, '\\d+'))

report$salmon1.pca.plot2 <- 
  make_square(ggplot(salmon1.pca.plotData2, aes(x=x, y=y, fill = time, shape = genotype, label = subject)) +
              theme_bw() +
              geom_point(size = 4, stroke = 0.5, color = 'black') +
              scale_shape_manual(name = 'Genotype', values = 21:25) +
              scale_fill_gradient(name = 'Time (days)', low = "white", high = "black") +
              labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
                   y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)')) +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))

ggsave(report$salmon1.pca.plot2, file = 'figures_and_tables/RNAseq1_no_Y664F.pdf')





# Second SALMON run PCA.
#--------------------------------------------------------------------------------------------------

# Transform the DESeq dataset into fold change values for each transcript.
salmon2.rld     <- rlog(salmon2.ddsTxi)    # Transform count data with regularized logarithm which returns log2 data normalized to library size.
salmon2.rlogMat <- assay(salmon2.rld)  # assay() extracts the matrix of normalized values.


# Perform principle component analysis.
i <- unname(apply(salmon2.rlogMat, 1, function(x) all(x < 15)))
salmon2.pca <- prcomp(t(salmon2.rlogMat[i,]), scale = FALSE, center = TRUE)


# Use the 'genotype_timePoint_donor' formatted data points to extract data for plots.
salmon2.pca.plotData          <- data.frame(s = row.names(salmon2.pca$x), x = salmon2.pca$x[,1], y = salmon2.pca$x[,2], z = salmon2.pca$x[,3])
salmon2.pca.plotData$genotype <- do.call(rbind, strsplit(salmon2.pca.plotData$s, '_'))[,1]
salmon2.pca.plotData$day      <- do.call(rbind, strsplit(salmon2.pca.plotData$s, '_'))[,2]
salmon2.pca.plotData$subject  <- do.call(rbind, strsplit(salmon2.pca.plotData$s, '_'))[,3]
salmon2.pca.plotData$subject  <- gsub('don', '', salmon2.pca.plotData$subject)
salmon2.pca.plotData          <- salmon2.pca.plotData[mixedorder(as.character(salmon2.pca.plotData$day)),]
salmon2.pca.plotData$day      <- factor(salmon2.pca.plotData$day, levels = unique(salmon2.pca.plotData$day))
salmon2.pca.plotData$genotype <- factor(salmon2.pca.plotData$genotype, levels = c('WT', 'KD', 'TrpM'))

report$salmon2.pca.plot <- 
  ggplot(salmon2.pca.plotData, aes(x=x, y=y, fill = genotype, shape = day, label = subject)) +
  theme_bw() +
  geom_point(size = 4, stroke = 0.75, color = 'black') +
  # geom_text(size = 5, nudge_x = 3, nudge_y = 3, show.legend = FALSE) +
  scale_shape_manual(name = 'Time point', values = 21:22) +
  scale_fill_manual(name = 'Transgene', values=c('blue', 'red', 'gray50')) +
  labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
       y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

ggsave(report$salmon2.pca.plot, file = 'figures_and_tables/RNAseq2.pdf')

save(list = ls(all.names = TRUE), file = 'savePoints/sp1.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
#-~-~-~-~o~-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~




# Create contrasts which contain the log2 fold change values and adjusted pvalues.
#--------------------------------------------------------------------------------------------------
RNAseq1_Y664F_D9_vs_WT_D9       <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d9", "WT_d9")))
RNAseq1_Y664F_D12_vs_WT_D12     <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d12", "WT_d12")))
RNAseq1_Y664F_D9_vs_KD_D9       <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d9", "KD_d9")))
RNAseq1_Y664F_D12_vs_KD_D12     <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d12", "KD_d12")))
RNAseq1_WT_D9_vs_WT_D6          <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d9", "WT_d6")))
RNAseq1_WT_D12_vs_WT_D6         <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d12", "WT_d6")))
RNAseq1_WT_D33_vs_WT_D6         <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d33", "WT_d6")))
RNAseq1_WT_D63_vs_WT_D6         <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d63", "WT_d6")))
RNAseq1_WT_D6_vs_KD_D6          <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d6", "KD_d6")))
RNAseq1_WT_D9_vs_KD_D9          <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d9", "KD_d9")))
RNAseq1_WT_D12_vs_KD_D12        <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d12", "KD_d12")))
RNAseq1_WT_D33_vs_KD_D12        <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d33", "KD_d12")))
RNAseq1_WT_D63_vs_KD_D12        <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d63", "KD_d12")))
RNAseq2_WT_D6_vs_KD_D6          <- addGeneNamesToContrast(results(salmon2.ddsTxi, contrast=c("repGrp","WT_D6","KD_D6")))
RNAseq2_TrpM_D6_vs_KD_D6        <- addGeneNamesToContrast(results(salmon2.ddsTxi, contrast=c("repGrp","TrpM_D6","KD_D6")))
RNAseq2_TrpM_D6_vs_WT_D6        <- addGeneNamesToContrast(results(salmon2.ddsTxi, contrast=c("repGrp","TrpM_D6","WT_D6")))
RNAseq2_WT_D9_vs_KD_D9          <- addGeneNamesToContrast(results(salmon2.ddsTxi, contrast=c("repGrp","WT_D9","KD_D9")))
RNAseq2_TrpM_D9_vs_KD_D9        <- addGeneNamesToContrast(results(salmon2.ddsTxi, contrast=c("repGrp","TrpM_D9","KD_D9")))
RNAseq2_TrpM_D9_vs_WT_D9        <- addGeneNamesToContrast(results(salmon2.ddsTxi, contrast=c("repGrp","TrpM_D9","WT_D9")))


# Build contrast tables
#------------------------------------------------------------------------------------------------00

RNAseq1_WT_vs_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_WT_D6_vs_KD_D6),   exp = 'WT D6',  timePoint = 'D6'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D9_vs_KD_D9),   exp = 'WT D9',  timePoint = 'D9'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D12_vs_KD_D12), exp = 'WT D12', timePoint = 'D12'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D33_vs_KD_D12), exp = 'WT D33', timePoint = 'D33'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D63_vs_KD_D12), exp = 'WT D63', timePoint = 'D63')) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = toupper(gene))


RNAseq1_WT_vs_d6_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_WT_D9_vs_WT_D6),  exp = 'D9',  timePoint = 'D9'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D12_vs_WT_D6), exp = 'D12', timePoint = 'D12'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D33_vs_WT_D6), exp = 'D33', timePoint = 'D33'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D63_vs_WT_D6), exp = 'D63', timePoint = 'D63')) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = toupper(gene))


RNAseq1_WT_vs_d12_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_WT_D12_vs_KD_D12), exp = 'D12', timePoint = 'D12'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D33_vs_KD_D12), exp = 'D33', timePoint = 'D33'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D63_vs_KD_D12), exp = 'D63', timePoint = 'D63')) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = toupper(gene))


RNAseq1_Y664F_vs_WT <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_Y664F_D9_vs_WT_D9),   exp = 'Y664F D9', timePoint = 'D9'),
                   dplyr::mutate(data.frame(RNAseq1_Y664F_D12_vs_WT_D12), exp = 'Y664F D12', timePoint = 'D12')) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = toupper(gene))



RNAseq1_Y664F_vs_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_Y664F_D9_vs_KD_D9),   exp = 'Y664F D9', timePoint = 'D9'),
                   dplyr::mutate(data.frame(RNAseq1_Y664F_D12_vs_KD_D12), exp = 'Y664F D12', timePoint = 'D12')) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = toupper(gene))



RNAseq2_WT_TrpM_vs_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq2_WT_D6_vs_KD_D6),   exp = 'WT D6',   timePoint = 'D6 WT'),
                   dplyr::mutate(data.frame(RNAseq2_TrpM_D6_vs_KD_D6), exp = 'TrpM D6', timePoint = 'D6 TrpM'),
                   dplyr::mutate(data.frame(RNAseq2_WT_D9_vs_KD_D9),   exp = 'WT D9',   timePoint = 'D9 WT'),
                   dplyr::mutate(data.frame(RNAseq2_TrpM_D9_vs_KD_D9), exp = 'TrpM D9', timePoint = 'D9 TrpM')) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(gene = toupper(gene))


RNAseq2_TrpM_vs_WT <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq2_TrpM_D6_vs_WT_D6), exp = 'D6', timePoint = 'D6'),
                   dplyr::mutate(data.frame(RNAseq2_TrpM_D9_vs_WT_D9), exp = 'D9', timePoint = 'D9')) %>%
  dplyr::mutate(exp = factor(exp, levels = c('D6', 'D9'))) %>%
  dplyr::filter(abs(log2FoldChange) <= 15) %>%
  tidyr::replace_na(list(padj = 1)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(gene = toupper(gene))



# Heatmaps
#--------------------------------------------------------------------------------------------------

ALK_ALCL_a    <- c('IL1RAP', 'BATF3', 'FBN1', 'TMEM158', 'IMPA2', 'FUT7', 'IL2RA', 'TMEM260', 'MYO10', 'MCAM', 'SLC12A8', 'NQO1', 'AGT', 'RHOC', 'PRF1', 'TEAD4', 'SPP1', 'CLEC3B', 'RUBCNL', 'G0S2', 'S100A9', 'TNFRSF8', 'PTPRG')
ALK_ALCL_b    <- c('IL10', 'TNFRSF8', 'GZMB', 'CD274', 'IL17A', 'IL22', 'JUNB', 'PRF1', 'CDC25A', 'LCK', 'BCL11B','LAT', 'LCP2', 'IL2RG', 'ZAP70', 'CD3E')
Tcell_TFs     <- c('TCF7', 'GATA3', 'LEF1', 'BCL11B', 'BACH2', 'NFATC3', 'HDAC7', 'SATB1', 'KLF3', 'ETS1', 'CAMK4', 'KLF7', 'MYC', 'THEMIS')
Tcell_sig     <- c('ZAP70', 'LCK', 'ITK', 'FYB1', 'NCK2', 'GRAP2', 'RASGRP1', 'FOXO1', 'PTPN6', 'LTB', 'MAL', 'ITPKB', 'LAT')
Tcell_recp    <- c('TRAC', 'TRBC2', 'TRBC1', 'TRGC1', 'CD3E', 'CD3D', 'CD3G', 'CD28', 'CD27', 'IL2RG', 'CCR7', 'CTLA4', 'CD2')
NonTcell_spec <- c('HHEX', 'JUNB', 'MEIS1', 'LAT2', 'LYN', 'NFE2', 'PAX5', 'MPO', 'GADD45A', 'MAFB', 'MEF2C', 'CD34', 'FLT3', 'GATA2')
EmbStemSpec   <- c('SOX2', 'SOX4', 'SOX5', 'TEAD4', 'ZIC2', 'ZIC5', 'ZEB2', 'DTX1', 'TWIST1', 'EDN1', 'PDGFA', 'CCND1', 'ETS2', 'EZH2', 'ZIC1') 
upInWnt       <- c('HHEX','JUNB', 'MEIS1','LAT2', 'LYN', 'NFE2', 'PAX5', 'MPO', 'GADD45A','MAFB', 'MEF2C', 'CD34', 'FLT3', 'GATA2', 'SOX2', 'SOX4', 'SOX5', 'TEAD4', 'ZIC2', 'ZIC5','ZEB2', 'DTX1', 'TWIST1','EDN1', 'PDGFA', 'CCND1', 'ETS2', 'EZH2', 'ZIC1')


# WT63vsKD12
Thymocyte_CD34aPLUS_CD1aMINUS <- c('LAT2','LYN','HHEX','TIGIT','LEF1','LCK','FYN','CAMK4','BCL11B')
Thymocyte_CD34aPLUS_CD1aPLUS  <- c('LAT2','LYN','RUNX2','HHEX','ITK','ZFP90','CD27','CTLA4','CD28')
Thymocyte_Pre_TCR             <- c('AURKA','MYB','EZH2','ETS1','TIGIT','RCAN3','CD28','ATM','ITGAL')
Thymocyte_DP_TCRa_b           <- c('EZH2','MYB','ATM','BACH2','RGS1','CD27','CTLA4','SERPINE2','KLF2')

report$Thymocyte_CD34aPLUS_CD1aMINUS_heatmap <- createGeneListHeatMap(subset(RNAseq1_WT_vs_KD, timePoint == 'D63'), Thymocyte_CD34aPLUS_CD1aMINUS, 9, 'figures_and_tables/RNAseq1_WT_vs_KD_D63_Thymocyte_CD34aPLUS_CD1aMINUS.pdf')
report$Thymocyte_CD34aPLUS_CD1aPLUS <- createGeneListHeatMap(subset(RNAseq1_WT_vs_KD, timePoint == 'D63'), Thymocyte_CD34aPLUS_CD1aPLUS, 9, 'figures_and_tables/RNAseq1_WT_vs_KD_D63_Thymocyte_CD34aPLUS_CD1aPLUS.pdf')
report$Thymocyte_Pre_TCR <- createGeneListHeatMap(subset(RNAseq1_WT_vs_KD, timePoint == 'D63'), Thymocyte_Pre_TCR, 9, 'figures_and_tables/RNAseq1_WT_vs_KD_D63_Thymocyte_Pre_TCR.pdf')
report$Thymocyte_DP_TCRa_b <- createGeneListHeatMap(subset(RNAseq1_WT_vs_KD, timePoint == 'D63'), Thymocyte_DP_TCRa_b, 9, 'figures_and_tables/RNAseq1_WT_vs_KD_D63_Thymocyte_DP_TCRa_b.pdf')

# Create heat maps for gene of interest sets where select sets share common scales.
report$ALK_ALCL_heatmap_a    <- createGeneListHeatMap(RNAseq1_WT_vs_KD, ALK_ALCL_a,    9,  'figures_and_tables/RNAseq1_WT_vs_KD_ALK_ALCL_heatmap_a.pdf')
report$ALK_ALCL_heatmap_b    <- createGeneListHeatMap(RNAseq1_WT_vs_KD, ALK_ALCL_b,    9,  'figures_and_tables/RNAseq1_WT_vs_KD_ALK_ALCL_heatmap_b.pdf')
report$Tcell_TFs_heatmap     <- createGeneListHeatMap(RNAseq1_WT_vs_KD, Tcell_TFs,     11, 'figures_and_tables/RNAseq1_WT_vs_KD_Tcell_TFs_heatmap.pdf')
report$Tcell_sig_heatmap     <- createGeneListHeatMap(RNAseq1_WT_vs_KD, Tcell_sig,     11, 'figures_and_tables/RNAseq1_WT_vs_KD_Tcell_sig_heatmap.pdf')
report$Tcell_recp_heatmap    <- createGeneListHeatMap(RNAseq1_WT_vs_KD, Tcell_recp,    11, 'figures_and_tables/RNAseq1_WT_vs_KD_Tcell_recp_heatmap.pdf')
report$NonTcell_spec_heatmap <- createGeneListHeatMap(RNAseq1_WT_vs_KD, NonTcell_spec, 11, 'figures_and_tables/RNAseq1_WT_vs_KD_NonTcell_spec_heatmap.pdf')
report$EmbStemSpec_heatmap   <- createGeneListHeatMap(RNAseq1_WT_vs_KD, EmbStemSpec,   11, 'figures_and_tables/RNAseq1_WT_vs_KD_EmbStemSpec_heatmap.pdf')


#------------------
RNAseq2_TrpM_vs_WT
RNAseq1_Y664F_vs_WT

a <- RNAseq2_TrpM_vs_WT[grep('D9', RNAseq2_TrpM_vs_WT$timePoint),]
b <- subset(RNAseq1_Y664F_vs_WT, timePoint == 'D9')
a$timePoint <- 'TrpM'
b$timePoint <- 'Y664F'
o <- bind_rows(a, b)
o$timePoint <- factor(o$timePoint, levels = c('Y664F', 'TrpM'))
o <- subset(o, gene %in% c(Tcell_sig, Tcell_TFs, Tcell_recp, upInWnt))
createGeneListHeatMap(o, Tcell_sig, ceiling(max(abs(o$log2FoldChange))),  'figures_and_tables/D9_WT_Y664F_TrpM_vs_KD_Tcell_sig_heatmap.pdf')

#--------------


# Create day 9 plots for WT, Y664F, and TrpM for select gene sets.
a <- RNAseq2_WT_TrpM_vs_KD[grep('D9', RNAseq2_WT_TrpM_vs_KD$timePoint),]
b <- subset(RNAseq1_Y664F_vs_KD, timePoint == 'D9')
b$timePoint <- 'D9 Y664F'
o <- bind_rows(a, b)
o$timePoint <- factor(sub('D9\\s+', '', o$timePoint), levels = c('WT', 'Y664F', 'TrpM'))
o <- subset(o, gene %in% c(Tcell_sig, Tcell_TFs, Tcell_recp, upInWnt))

report$D9_WT_Y664F_TrpM_vs_KD_Tcell_sig_heatmap  <- createGeneListHeatMap(o, Tcell_sig, ceiling(max(abs(o$log2FoldChange))),  'figures_and_tables/D9_WT_Y664F_TrpM_vs_KD_Tcell_sig_heatmap.pdf')
report$D9_WT_Y664F_TrpM_vs_KD_Tcell_TFs_heatmap  <- createGeneListHeatMap(o, Tcell_TFs, ceiling(max(abs(o$log2FoldChange))),  'figures_and_tables/D9_WT_Y664F_TrpM_vs_KD_Tcell_TFs_heatmap.pdf')
report$D9_WT_Y664F_TrpM_vs_KD_Tcell_recp_heatmap <- createGeneListHeatMap(o, Tcell_recp, ceiling(max(abs(o$log2FoldChange))), 'figures_and_tables/D9_WT_Y664F_TrpM_vs_KD_Tcell_recp_heatmap.pdf')
report$D9_WT_Y664F_TrpM_vs_KD_upInWNT_heatmap    <- createGeneListHeatMap(o, upInWnt, ceiling(max(abs(o$log2FoldChange))),    'figures_and_tables/D9_WT_Y664F_TrpM_vs_KD_upInWNT_heatmap.pdf')



# RNAseq2  WT / TrpM vs KD
# Determine which genes should be used - make sure that selected genes are seen across all time points.
genes <- dplyr::group_by(RNAseq2_WT_TrpM_vs_KD, gene) %>%
         dplyr::mutate(exps = n_distinct(exp), colsum = sum(log2FoldChange)) %>%
         dplyr::ungroup() %>%
         dplyr::filter(! is.na(padj) & abs(log2FoldChange) >= 3 & exps == n_distinct(exp) & padj <= 1e-70) %>%
         dplyr::arrange(desc(colsum)) %>%
         dplyr::pull(gene) %>%
         unique()

report$RNAseq2_WT_TrpM_vs_KD.plot1 <- createGeneListHeatMap(RNAseq2_WT_TrpM_vs_KD, genes[1:38],  11, 'figures_and_tables/RNAseq2_WT_TrpM_vs_KD.plot_part1.pdf', orderByFoldChange = FALSE)
report$RNAseq2_WT_TrpM_vs_KD.plot2 <- createGeneListHeatMap(RNAseq2_WT_TrpM_vs_KD, genes[39:75], 11, 'figures_and_tables/RNAseq2_WT_TrpM_vs_KD.plot_part2.pdf', orderByFoldChange = FALSE)




# RNAseq2 SWI/SNF
#--------------------------------------------------------------------------------------------------
genes    <- unique(RNAseq2_WT_TrpM_vs_KD$gene[grepl('^SMARC|^ARID1|^ACTL6', RNAseq2_WT_TrpM_vs_KD$gene)])
g        <- subset(RNAseq2_WT_TrpM_vs_KD, gene %in% genes)
report$SWI_SNF_genes_RNAseq <- createGeneListHeatMap(g, unique(g$gene), 6, 'figures_and_tables/RNAseq2_WT_TrpM_vs_KD_SWI_SNF.pdf')



# Compare WT transcription to earliest WT profile.
#--------------------------------------------------------------------------------------------------
RNAseq1_WT_vs_d6_KD.genes <- subset(data.frame(RNAseq1_WT_D63_vs_WT_D6), padj <= 1e-25 & abs(log2FoldChange) >= 3)$gene
report$RNAseq1.pval.WT_vs_earlyWT <- createGeneListHeatMap(RNAseq1_WT_vs_d6_KD, RNAseq1_WT_vs_d6_KD.genes, 8, 'figures_and_tables/RNAseq1_WT_vs_d6_KD.pdf')



# RNAseq1: Compare Y664F to WT. 
#--------------------------------------------------------------------------------------------------
RNAseq1_Y664F_vs_WT.genes <- 
  dplyr::group_by(RNAseq1_Y664F_vs_WT, gene) %>%
  dplyr::mutate(exps = n_distinct(as.character(exp))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(abs(log2FoldChange) >= 2 & exps == n_distinct(exp) & padj <= 1e-03) %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  unique()

report$RNAseq1.Y664F_vs_WT <- createGeneListHeatMap(RNAseq1_Y664F_vs_WT, RNAseq1_Y664F_vs_WT.genes, 9, 'figures_and_tables/RNAseq1.Y664F_vs_WT.pdf')




# RNAseq2: Compare TrpM to WT. 
#--------------------------------------------------------------------------------------------------
RNAseq2_TrpM_vs_WT.genes <- 
  dplyr::group_by(RNAseq2_TrpM_vs_WT, gene) %>%
  dplyr::mutate(exps = n_distinct(exp)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(! is.na(padj) & abs(log2FoldChange) >= 4 & exps == n_distinct(exp) & padj <= 1e-15) %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  unique()

report$RNAseq2.TrpM_vs_WT <- createGeneListHeatMap(RNAseq2_TrpM_vs_WT, RNAseq2_TrpM_vs_WT.genes, 8, 'figures_and_tables/RNAseq2.TrpM_vs_WT.pdf')


# Volcano plots
#--------------------------------------------------------------------------------------------------
report$RNAseq1_WT_D6_vs_KD_D6_volcano   <- createVolcanoPlot(RNAseq1_WT_D6_vs_KD_D6,   'WT D6 vs KD D6',   'figures_and_tables/RNAseq1_WT_D6_vs_KD_D6_volcano.pdf')
report$RNAseq1_WT_D9_vs_KD_D9_volcano   <- createVolcanoPlot(RNAseq1_WT_D9_vs_KD_D9,   'WT D9 vs KD D9',   'figures_and_tables/RNAseq1_WT_D9_vs_KD_D9_volcano.pdf')
report$RNAseq1_WT_D12_vs_KD_D12_volcano <- createVolcanoPlot(RNAseq1_WT_D12_vs_KD_D12, 'WT D12 vs KD D12', 'figures_and_tables/RNAseq1_WT_D12_vs_KD_D12_volcano.pdf')
report$RNAseq1_WT_D33_vs_KD_D12_volcano <- createVolcanoPlot(RNAseq1_WT_D33_vs_KD_D12, 'WT D33 vs KD D12', 'figures_and_tables/RNAseq1_WT_D33_vs_KD_D12_volcano.pdf')
report$RNAseq1_WT_D63_vs_KD_D12_volcano <- createVolcanoPlot(RNAseq1_WT_D63_vs_KD_D12, 'WT D63 vs KD D12', 'figures_and_tables/RNAseq1_WT_D63_vs_KD_D12_volcano.pdf')



save(list = ls(all.names = TRUE), file = 'savePoints/sp2.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
#-~-~-~-~o~-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~




# Add entrez and STRINGdb ids to the RNAseq data.
#--------------------------------------------------------------------------------------------------
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "entrezgene_accession"), mart = ensembl)
mapping$entrezgene <- mapping$entrezgene_id

# The mapping object is saved in 'data/bioMartMapping.RData' in case biomart goes off-line.

RNAseq1_WT_vs_KD      <- addAndExtendSTRINGids(RNAseq1_WT_vs_KD)
RNAseq1_WT_vs_d6_KD   <- addAndExtendSTRINGids(RNAseq1_WT_vs_d6_KD)
RNAseq1_WT_vs_d12_KD  <- addAndExtendSTRINGids(RNAseq1_WT_vs_d12_KD)
RNAseq2_WT_TrpM_vs_KD <- addAndExtendSTRINGids(RNAseq2_WT_TrpM_vs_KD)


# GSEA
keggGenes <- gage::kegg.gsets(species = "hsa", id.type = "kegg")

o <- subset(RNAseq1_WT_vs_KD, timePoint == 'D9')
ranks <- o$log2FoldChange
names(ranks) <- o$entrezgene
barplot(sort(ranks, decreasing = T))

ranks <- ranks[! is.na(names(ranks))]
fgseaRes <- fgsea(keggGenes$kg.sets, ranks, minSize = 15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=10)

barplot(sort(ranks, decreasing = T))
plotEnrichment(keggGenes$kg.sets[["hsa04062 Chemokine signaling pathway"]], ranks)
plotEnrichment(keggGenes$kg.sets[["hsa04630 Jak-STAT signaling pathway"]], ranks)






# KEGG pathway enrichments
#--------------------------------------------------------------------------------------------------
# Detach the RMySQL package because it interferes with STRINGdb
detach('package:RMySQL', unload = TRUE, character.only = TRUE)


invisible(lapply(split(RNAseq1_WT_vs_KD, RNAseq1_WT_vs_KD$timePoint), function(x){
  o <- list()
  o$title <- paste(x$timePoint[1], 'WT vs', 
                 ifelse(as.integer(str_extract(x$timePoint[1], '\\d+')) > 12, 'D12 KD', paste0(x$timePoint[1], ' KD')))
  
  file <- paste0('KEGG_pathways_', gsub(' ', '_', o$title))

  o$table     <- createKEGGenrichmentTable(subset(x, ! is.na(log2FoldChange) & ! is.na(STRING_id) & padj <= 1e-3 & abs(log2FoldChange) >= 2))
  o$plot      <- termEnrichmentPlot(o$table, 10, o$title, file = paste0('figures_and_tables/RNAseq1_', file, '.pdf'))
  o$up.plot   <- termEnrichmentPlot(o$table, 10, o$title, file = paste0('figures_and_tables/RNAseq1_', file, '.up.pdf'), dir = 'up')
  o$down.plot <- termEnrichmentPlot(o$table, 10, o$title, file = paste0('figures_and_tables/RNAseq1_', file, '.down.pdf'), dir = 'down')
  o
}))



# GO term enrichment
#--------------------------------------------------------------------------------------------------
report$WT_D9_vs_KD_D9.GO <- 
  subset(string_db$get_enrichment(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'WT D9' & 
                                    ! is.na(log2FoldChange) &  
                                    ! is.na(STRING_id) & 
                                    padj <= 1e-3 & 
                                    abs(log2FoldChange) >= 2)$STRING_id, category = 'Process', methodMT = "fdr", iea = FALSE),
         pvalue_fdr <= 0.05)

write.xlsx(report$WT_D9_vs_KD_D9.GO, file = file.path('figures_and_tables', 'WT_D9_vs_KD_D9.GO_terms.xlsx'), col.names = TRUE, row.names = FALSE)

report$TrpM_D9_vs_KD_D9.GO <- 
  subset(string_db$get_enrichment(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'TrpM D9' & 
                                    ! is.na(log2FoldChange) &  
                                    ! is.na(STRING_id) & 
                                    padj <= 1e-3 & 
                                    abs(log2FoldChange) >= 2)$STRING_id, category = 'Process', methodMT = "fdr", iea = FALSE), 
         pvalue_fdr <= 0.05)


write.xlsx(report$TrpM_D9_vs_KD_D9.GO, file = file.path('figures_and_tables', 'TrpM_D9_vs_KD_D9.GO_terms.xlsx'), col.names = TRUE, row.names = FALSE)

save(list = ls(all.names = TRUE), file = 'savePoints/sp3.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)
#-~-~-~-~o~-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~o-~-~-~-~



# Integration site data
#-------------------------------------------------------------------------------------------------------

report$intSites <- intSites; rm(intSites)
report$intSites$genotype <- ifelse(grepl('Y664F', report$intSites$patient), 'Y664F', 'WT')


# Create a list of relative abundance plots for each subject.

report$sampleAbundancePlotsData_n <- 50
d <- data.frame(report$intSites)
o <- split(d, d$patient)
report$sampleAbundancePlotsData <- o[order(sapply(o, function(x){ t <- sort(x$timePointDays); t[1] } ))]

report$sampleAbundancePlots <- lapply(report$sampleAbundancePlotsData, function(x){
  x <- x[rev(order(x$relAbund)),]
  x$geneLabel <- paste0(x$nearestFeature, '\n', x$posid)
  if(length(unique(x$posid)) < report$sampleAbundancePlotsData_n) report$sampleAbundancePlotsData_n<- length(unique(x$posid))
  sites <- x$geneLabel[1:report$sampleAbundancePlotsData_n]
  
  df <- do.call(rbind, lapply(split(x, x$timePoint), function(x2){
    x2 <- subset(x2, geneLabel %in% sites)
    x2 <- x2[rev(order(x2$relAbund)),]
    data.frame(patient=x2$patient[1],
               timePoint=x2$timePoint[1],
               relAbund=c(100-sum(x2$relAbund), x2$relAbund),
               sites=c('LowAbund', x2$geneLabel),
               stringsAsFactors = FALSE)
  }))
  
  df$sites     <- factor(df$sites, levels=unique(df$sites))
  df$timePoint <- factor(df$timePoint, levels=mixedsort(unique(df$timePoint)))
  
  ggplot() + 
    theme_bw() +
    geom_bar(data=df, aes(timePoint, relAbund/100, fill=sites), stat = "identity") + 
    guides(fill=FALSE) +
    labs(x='', y='') +
    ggtitle(x$patient[1]) +
    scale_fill_manual(values=c('#DCDCDC', colorRampPalette(brewer.pal(12, "Paired"))(report$sampleAbundancePlotsData_n))) +
    theme(plot.title = element_text(size = 12)) +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    scale_y_continuous(labels = scales::percent) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))
})


write.table(bind_rows(report$sampleAbundancePlotsData) %>% 
            dplyr::select(seqnames, start, strand, patient, cellType, timePoint, estAbund, relAbund, nearestFeature, genotype), 
            file = file.path('figures_and_tables', 'intSites_relativeAbund_plotData.csv'), sep = ',', col.names = TRUE, row.names = FALSE)


# Create list of patients with multiple time points.
report$subjectsWithMultTimePoints <- 
  dplyr::group_by(data.frame(report$intSites), patient) %>% 
  dplyr::summarise(tps = n_distinct(timePoint)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(tps > 1) %>% 
  dplyr::select(patient) %>% 
  unlist() %>% unname()



# Map STRINGdb ids to intSite nearest genes.
d <- string_db$map(data.frame(report$intSites), "nearestFeature", removeUnmappedRows = FALSE)
i <- which(is.na(d$STRING_id));  message(paste0('gene ids not mapped to STRINGdb: ', round((length(i) / nrow(d))*100, 2), '%'))
d[i,]$STRING_id <- report$string_db.alt.aliases[match(toupper(d[i,]$nearestFeature), toupper(report$string_db.alt.aliases$ids)),]$STRING_id
i <- which(is.na(d$STRING_id));  message(paste0('gene ids not mapped to STRINGdb: ', round((length(i) / nrow(d))*100, 2), '%'))
report$intSites$STRING_id <- d$STRING_id



# Create an integration frequency bivariate plot for WT subjects.
report$preVspostFreqplot.WT <- preVsPostFreqPlot(subset(report$intSites, patient %in% c('pNA92', 'pNA93', 'pNA85')),
                                                 14, 
                                                 distCutOff = 50000,
                                                 readRDS('data/humanOncoGenesList.rds'), 
                                                 nGenesToLabel = 5,
                                                 mustLabel = c('ARID1A'))


ggsave(report$preVspostFreqplot.WT, file = 'figures_and_tables/preVspostFreqplot.WT.pdf')


# Create an integration frequency bivariate plot for Y664F subjects.
report$preVspostFreqplot.Y664F <- preVsPostFreqPlot(subset(report$intSites, patient %in% c('pNA92_Y664F', 'pNA93_Y664F', 'pNA85_Y664F')),
                                                    14, 
                                                    readRDS('data/humanOncoGenesList.rds'), 
                                                    nGenesToLabel = 5,
                                                    mustLabel = c('ARID1A'))


ggsave(report$preVspostFreqplot.Y664F, file = 'figures_and_tables/preVspostFreqplot.Y664F.pdf')


# Create UCSD tracks.
report$intsites.WT       <- subset(report$intSites, ! grepl('Y664F',patient))
report$intsites.Y664F    <- subset(report$intSites, grepl('Y664F',patient))
report$intsites.LS       <- subset(report$intSites, patient %in% report$subjectsWithMultTimePoints)
report$intsites.LS_WT    <- subset(report$intSites, ! grepl('Y664F',patient) & patient %in% report$subjectsWithMultTimePoints) 
report$intsites.LS_Y664F <- subset(report$intSites, grepl('Y664F',patient) & patient %in% report$subjectsWithMultTimePoints)


createUCSCintSiteAbundTrack(report$intSites$posid, report$intSites$estAbund, report$intSites$patient, title = 'All_samples', visbility = 1, outputFile = 'allSamples.ucsd')
createUCSCintSiteAbundTrack(report$intsites.WT$posid, report$intsites.WT$estAbund, report$intsites.WT$patient, title = 'All_WT_samples', visbility = 0, outputFile = 'allWTsamples.ucsd')
createUCSCintSiteAbundTrack(report$intsites.Y664F$posid, report$intsites.Y664F$estAbund, report$intsites.Y664F$patient, title = 'All_Y664F_samples', visbility = 0, outputFile = 'allY664Fsamples.ucsd')
system(paste('cat', paste(list.files(pattern = '*.ucsd'), collapse = ' '), '> NPM_ALK.ucsd'))
system(paste0('scp NPM_ALK.ucsd  microb120:/usr/share/nginx/html/UCSC/everett/NPM_ALK/'))
invisible(file.remove(list.files(pattern = '\\.ucsd$')))



d <- data.frame(report$intSites)
d <- subset(d, patient %in% report$subjectsWithMultTimePoints)

a <- unique(subset(d, timePointDays <= 14)$posid)
b <- unique(subset(d, timePointDays > 14)$posid)
report$type1_posids <- a[! a %in% b]

report$percentClonesOnlySeenAtStart <- sprintf("%.2f%%", (n_distinct(report$type1_posids) / n_distinct(d$posid))*100)

d$posid <- paste(d$patient, d$posid)

# All sites need a start point.
d <- bind_rows(lapply(split(d, d$patient), function(x){
     all_posids <- unique(x$posid)
     startPosids <- subset(x, timePoint == 'D12')$posid
     
     # There is a D12 time point and all intSites are not found in it.
     if(length(startPosids) > 0 & ! all(all_posids %in% startPosids)){
         missing_posids <- all_posids[which(! all_posids %in% startPosids)]
         o <- subset(x, posid %in% missing_posids)
     } else {
       o <- x
     }
     
     o <- o[! duplicated(o$posid),]
     o$estAbund <- 0
     o$timePoint <- 'D12'
     o$timePointDays <- 12
     o$timePointMonths <- 0.395
     bind_rows(x, o)
}))

# Add zero points for longitudinal trials.
d <- bind_rows(lapply(split(d, d$patient), function(x){
  all_posids <- unique(x$posid)
  bind_rows(lapply(unique(x$timePoint), function(tp){
    x2 <- subset(x, timePoint == tp)
    if(! all(all_posids %in% x2$posid)){
      missing_posids <- all_posids[which(! all_posids %in% x2$posid)]
      o <- subset(x, posid %in% missing_posids)
      o <- o[! duplicated(o$posid),]
      o$estAbund <- 0
      o$timePoint <- tp
      o$timePointDays <- timePoint2numeric(tp, interval = 'days')
      o$timePointMonths <- timePoint2numeric(tp, interval = 'months')
      return(bind_rows(x2, o))
    } else {
      return(x2)
    }
  }))
}))


d.maxEstAbund <- 
  dplyr::group_by(d, patient, posid) %>%
  dplyr::arrange(timePointDays) %>%
  dplyr::mutate(lastTimePointGreatestAbund = ifelse(max(estAbund) == estAbund[n()], TRUE, FALSE)) %>%
  dplyr::select(patient, posid, estAbund, timePointDays, genotype, nearestFeature, nearestFeatureDist, STRING_id, lastTimePointGreatestAbund) %>%
  dplyr::top_n(1, wt = estAbund) %>%
  dplyr::ungroup() 



p <- subset(d, posid %in% subset(d.maxEstAbund, estAbund >= 25)$posid)
p$class <- 'Decreasing'
p[which(p$posid %in% subset(d.maxEstAbund, lastTimePointGreatestAbund == TRUE)$posid),]$class <- 'Increasing'


# Hot fix
p[which(p$posid == 'pNA92 chr2+29127458'),]$class <- 'Increasing'

report$clonAbundanceTrajectories.colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_distinct(p$posid))
report$clonAbundanceTrajectories <- 
  ggplot(p, aes(timePointDays, estAbund, color = posid, group = posid, shape = genotype)) +
  theme_bw() +
  geom_point(size = 2.5) +
  geom_line() +
  scale_color_manual(values =  report$clonAbundanceTrajectories.colors) +
  guides(color = FALSE) +
  labs(x = 'Days', y = 'Abundance') + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) +
  facet_grid(class~.)

report$type2_posids <- unique(subset(p, class == 'Decreasing')$posid)
report$type3_posids <- unique(subset(p, timePointDays > 50 & estAbund <= 250 & class == 'Increasing')$posid)
report$type4_posids <- unique(subset(p, timePointDays > 50 & estAbund > 250 & class == 'Increasing')$posid)

report$type2_posids <- unique(unlist(lapply(stringr::str_split(report$type2_posids, '\\s'), '[[', 2)))
report$type3_posids <- unique(unlist(lapply(stringr::str_split(report$type3_posids, '\\s'), '[[', 2)))
report$type4_posids <- unique(unlist(lapply(stringr::str_split(report$type4_posids, '\\s'), '[[', 2)))

# Remove patient ids from posids.
d$posid <- unlist(lapply(stringr::str_split(d$posid, '\\s'), '[[', 2))
p$posid <- unlist(lapply(stringr::str_split(p$posid, '\\s'), '[[', 2))

oncogenes <- toupper(c(readRDS('data/humanOncoGenesList.rds'), 'TET2,TET2-AS1'))

# Create a list of enriched KEGG terms from the D9 WT vs KD RNAseq experiment
k <- string_db$get_term_proteins(report$WT_D9_vs_KD_D9.KEGG$term_id)$STRING_id

type2_string_ids <- dplyr::group_by(subset(p, posid %in% report$type2_posids), posid) %>% dplyr::summarise(id = STRING_id[1]) %>% dplyr::ungroup() %>% dplyr::select(id) %>% unlist() %>% unname()
type3_string_ids <- dplyr::group_by(subset(p, posid %in% report$type3_posids), posid) %>% dplyr::summarise(id = STRING_id[1]) %>% dplyr::ungroup() %>% dplyr::select(id) %>% unlist() %>% unname()
type4_string_ids <- dplyr::group_by(subset(p, posid %in% report$type4_posids), posid) %>% dplyr::summarise(id = STRING_id[1]) %>% dplyr::ungroup() %>% dplyr::select(id) %>% unlist() %>% unname()

type2_genotypes <- dplyr::group_by(subset(p, posid %in% report$type2_posids), posid) %>% dplyr::summarise(genotype = genotype[1]) %>% dplyr::ungroup() %>% dplyr::select(genotype) %>% table() %>% paste(collapse = ' / ')
type3_genotypes <- dplyr::group_by(subset(p, posid %in% report$type3_posids), posid) %>% dplyr::summarise(genotype = genotype[1]) %>% dplyr::ungroup() %>% dplyr::select(genotype) %>% table() %>% paste(collapse = ' / ')
type4_genotypes <- dplyr::group_by(subset(p, posid %in% report$type4_posids), posid) %>% dplyr::summarise(genotype = genotype[1]) %>% dplyr::ungroup() %>% dplyr::select(genotype) %>% table() %>% paste(collapse = ' / ')

createCloneList <- function(posidList){
  subset(p, posid %in% posidList) %>% 
  dplyr::group_by(posid) %>% 
  top_n(1, wt = estAbund) %>%
  dplyr::select(posid, estAbund, nearestFeature, nearestFeatureDist) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(oncogene = toupper(nearestFeature) %in% oncogenes) %>% 
  dplyr::arrange(desc(estAbund))
}


p[which(p$posid == 'chr7+50229821'),]$nearestFeature <- '* AC020743.2'
p[which(p$posid == 'chr7+50229821'),]$nearestFeatureDist <- 0
p[which(p$posid == 'chr7+17293342'),]$nearestFeature <- '* anti-AHR'
p[which(p$posid == 'chr7+17293342'),]$nearestFeatureDist <- 0
p[which(p$posid == 'chr19+9648087'),]$nearestFeatureDist <- 0
p[which(p$posid == 'chr4-173365250'),]$nearestFeature <- '* AC097534.2'
p[which(p$posid == 'chr4-173365250'),]$nearestFeatureDist <- 0
p[which(p$posid == 'chr8+100404494'),]$nearestFeature <- '* AP003472.1'
p[which(p$posid == 'chr8+100404494'),]$nearestFeatureDist <- 0

report$type2_clones <- createCloneList(report$type2_posids)
report$type3_clones <- createCloneList(report$type3_posids)
report$type4_clones <- createCloneList(report$type4_posids)


report$clonAbundanceTrajectoriesTable <- data.frame(Class = c('Type 2', 'Type 3', 'Type 4'), 
                'Clones' = c(n_distinct(report$type2_posids),
                             n_distinct(report$type3_posids),
                             n_distinct(report$type4_posids)),
                '% Nearest gene in RNAseq enriched pathways' = c(sprintf("%.2f%%", (sum(type2_string_ids %in% k) / n_distinct(type2_string_ids))*100),
                                                                 sprintf("%.2f%%", (sum(type3_string_ids %in% k) / n_distinct(type3_string_ids))*100),
                                                                 sprintf("%.2f%%", (sum(type4_string_ids %in% k) / n_distinct(type4_string_ids))*100)), 
                'WT / Y664F' = c(type2_genotypes, type3_genotypes, type4_genotypes),
                check.names = FALSE)

report$type2.KEGG <- string_db$get_enrichment(type2_string_ids, category = 'KEGG', methodMT = "fdr", iea = FALSE)
report$type2.KEGG <- subset(report$type2.KEGG, pvalue_fdr <= 0.05)

report$type3.KEGG <- string_db$get_enrichment(type3_string_ids, category = 'KEGG', methodMT = "fdr", iea = FALSE)
report$type3.KEGG <- subset(report$type3.KEGG, pvalue_fdr <= 0.05)

report$type4.KEGG <- string_db$get_enrichment(type4_string_ids, category = 'KEGG', methodMT = "fdr", iea = FALSE)
report$type4.KEGG <- subset(report$type4.KEGG, pvalue_fdr <= 0.05)


report$intSites.maxEstAbund <- 
  dplyr::group_by(data.frame(report$intSites), patient, posid) %>%
  dplyr::arrange(timePointDays) %>%
  dplyr::mutate(lastTimePointGreatestAbund = ifelse(max(estAbund) == estAbund[n()], TRUE, FALSE)) %>%
  dplyr::select(patient, posid, estAbund, timePointDays, genotype, nearestFeature, nearestFeatureDist, STRING_id, lastTimePointGreatestAbund) %>%
  dplyr::top_n(1, wt = estAbund) %>%
  dplyr::ungroup() 


genesFoundAcrossMultSubject <-  dplyr::filter(data.frame(report$intSites.maxEstAbund), estAbund >= 25) %>%
  dplyr::group_by(nearestFeature) %>%
  dplyr::summarise(Y664F.subjects = sum(grepl('Y664F', unique(patient))),
                   WT.subjects = n_distinct(patient) - Y664F.subjects,
                   n = Y664F.subjects + WT.subjects) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n >= 3) %>%
  dplyr::arrange(desc(n)) 

report$genesFoundAcrossMultSubject <- bind_rows(lapply(split(genesFoundAcrossMultSubject, genesFoundAcrossMultSubject$nearestFeature), function(x){
  x$WT.maxAbundance <- max(subset(report$intSites.maxEstAbund, nearestFeature == x$nearestFeature & genotype == 'WT' & estAbund >= 25)$estAbund)
  if(is.infinite(x$WT.maxAbundance)) x$WT.maxAbundance <- NA
  x$Y664F.maxAbundance <- max(subset(report$intSites.maxEstAbund, nearestFeature == x$nearestFeature & genotype == 'Y664F' & estAbund >= 25)$estAbund)
  if(is.infinite(x$Y664F.maxAbundance)) x$Y664F.maxAbundance <- NA
  x
}))

report$genesFoundAcrossMultSubject <- 
  dplyr::select(report$genesFoundAcrossMultSubject, nearestFeature, WT.subjects, Y664F.subjects, WT.maxAbundance, Y664F.maxAbundance, n) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::select(-n)

saveRDS(report, file = 'report.rds')

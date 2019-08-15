library(tidyverse)
library(parallel)
library(RMySQL)
library(gt23)
library(parallel)
library(STRINGdb)
library(biomaRt)
library(stringr)
library(DESeq2)
library(gtools)
library(pathview)
library(xlsx)
source('./lib.R')
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE, java.parameters = "-Xmx4g") 

# Create a report list object to store data need for final report.
report <- list()
report$CPUs               <- 40
report$STRINGdb_dataFiles <- 'data/STRINGdb' # Created on the fly if not defined.
report$HGNC_dataFile      <- 'data/HGNC.txt'
report$salmonCommand      <- '/home/everett/software/salmon-0.11.3/bin/salmon quant -p8 -i /home/everett/data/sequenceDatabases/Salmon/GRCh38.gencode29/ -l A'
report$runSalmonRun1      <- FALSE
report$runSalmonRun2      <- FALSE

# Create HGNC aliases for better KEGG and GO enrichments.
#--------------------------------------------------------------------------------------------------

cluster <- makeCluster(report$CPUs)

# STRINGdb will be used for pathway and GO term enrichment.
# This approach is limited by the mapping of transcripts to genes and the mapping of genes to annotations.
# STRINGdb understands a number of gene aliases but the alias list is incompleted.
# Here we use the previous gene ids and gene synonyms provided by HGNC to expand STRINGdb's gene alias -> protein mappings.

HGNC <- read.table(report$HGNC_dataFile, sep = '\t', comment.char = '', quote = '', header = TRUE)

string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory = report$STRINGdb_dataFiles)


# For each vector of gene ids, determine which ones are in the current STRINGdb alias list.
# If any of the ids for a gene are in the alias list, assoicate all ids with the STRINGdb protein id
# but only if all of the aliases only return a single string ID.

if(! file.exists('data/string_db.alt.aliases.rds')){

  HGNC.ids <- 
    dplyr::rowwise(HGNC) %>%
    dplyr::summarise(ids = list(gsub('\\s', '', 
                                     toupper(c(Approved.Symbol, 
                                               unlist(strsplit(Previous.Symbols, ',')), 
                                               unlist(strsplit(Synonyms, ','))))))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n = ntile(1:n(), report$CPUs))
  
  STRINGdb_dataFiles <- report$STRINGdb_dataFiles
  clusterExport(cluster, varlist = c('STRINGdb_dataFiles', 'HGNC.ids'))

  report$string_db.alt.aliases <- bind_rows(parLapply(cluster, split(HGNC.ids, HGNC.ids$n), function(x){
    library(STRINGdb)
    library(dplyr)
    remove_NA <- function(x) x[!is.na(x)]
    
    string_db <- STRINGdb$new(version = "10", species = 9606, score_threshold = 0, input_directory = STRINGdb_dataFiles)
    string_db.aliases <- string_db$get_aliases()
    string_db.aliases$alias <- gsub('\\s', '', toupper(string_db.aliases$alias))
    
    dplyr::rowwise(x) %>%
    dplyr::mutate(stringIDS = list(unique(remove_NA(string_db.aliases[match(unlist(ids), string_db.aliases$alias),]$STRING_id)))) %>%
    dplyr::mutate(STRING_id = ifelse(length(stringIDS) == 1, unlist(stringIDS), NA)) %>%
    dplyr::ungroup() %>%
    dplyr::select(ids, STRING_id) %>% 
    dplyr::filter(! is.na(STRING_id)) %>% 
    tidyr::unnest(ids)
  }))

  saveRDS(report$string_db.alt.aliases, file = 'data/string_db.alt.aliases.rds')
} else {
  report$string_db.alt.aliases <- readRDS('data/string_db.alt.aliases.rds')
}


# Run SALMON data set1 if requested.
if(report$runSalmonRun1){
  if(! dir.exists('salmonOutput')) dir.create('salmonOutput')
  salmonRuns1 <- read.table('data/RNAseq_data1.tsv', header = TRUE, sep = '\t') %>%
                 dplyr::rowwise() %>% 
                 dplyr::select(filename, barcode) %>% 
                 dplyr::mutate(command = paste0(report$salmonCommand, 
                                         ' -1 ', list.files('data/RNAseq_data1', pattern = barcode, full.names = TRUE)[1], 
                                         ' -2 ', list.files('data/RNAseq_data1', pattern = barcode, full.names = TRUE)[2],
                                         ' -o salmonOutput/', filename)) %>%
                 dplyr::ungroup()
  invisible(sapply(salmonRuns1$command, system))
}


# Run SALMON data set2 if requested.
if(report$runSalmonRun2){
  if(! dir.exists('salmonOutput')) dir.create('salmonOutput')
  salmonRuns2 <- read.table('data/RNAseq_data2.tsv', header = TRUE, sep = '\t') %>%
                 dplyr::rowwise() %>% 
                 dplyr::mutate(command = paste0(report$salmonCommand, 
                                                ' -r ', dataPath,
                                                ' -o salmonOutput/', sample)) %>%
                 dplyr::ungroup()
  invisible(sapply(salmonRuns2$command, system))
}


# Import Salmon RNAseq results.
#--------------------------------------------------------------------------------------------------

# Import and name salmon output files.
# The two RNAseq runs, which used different strategies, can be separated by participant numbers.

files <- list.files(path = 'salmonOutput', pattern = 'quant.sf', recursive = TRUE, full.names = TRUE)
names(files) <- unlist(lapply(strsplit(files, '\\/'), '[[', 2))

salmon1.files <- files[as.integer(str_extract(names(files), '\\d+$')) <= 3]

# Exclude all day 9 results from subject 4 per Jan's instructions.
salmon2.files <- files[as.integer(str_extract(names(files), '\\d+$')) > 3]
salmon2.files <- salmon2.files[! grepl('D9_don4', names(salmon2.files))]


# Create RNAseq data import objects.
salmon1.txi <- importSalmon(salmon1.files)
salmon2.txi <- importSalmon(salmon2.files)


# Create DEseq2 data objects which can be used to retrieve counts and create contrasts.
salmon1.ddsTxi <- processSalmon(salmon1.txi)
salmon2.ddsTxi <- processSalmon(salmon2.txi)




# RNAseq PCA plots
#--------------------------------------------------------------------------------------------------

# First SALMON run PCA.


salmon1.rld <- rlog(salmon1.ddsTxi) # Transform count data with regularized logarithm which returns log2 data normalized to library size. 
salmon1.rlogMat <- assay(salmon1.rld) # assay() extracts the matrix of normalized values.

# i <- unname(apply(salmon1.rlogMat, 1, function(x) all(x < 14))) # Exclude outlier data points which were not caught by SALMON.
# salmon1.pca <- prcomp(t(salmon1.rlogMat[i,]), scale = FALSE, center = TRUE)
salmon1.pca <- prcomp(t(salmon1.rlogMat), scale = FALSE, center = TRUE)

# Use the 'genotype_timePoint_donor' formatted data points to extract data for the plots
salmon1.pca.plotData          <- data.frame(s = row.names(salmon1.pca$x), x = salmon1.pca$x[,1], y = salmon1.pca$x[,2], z = salmon1.pca$x[,3])
salmon1.pca.plotData$genotype <- do.call(rbind, strsplit(salmon1.pca.plotData$s, '_'))[,1]
salmon1.pca.plotData$day      <- do.call(rbind, strsplit(salmon1.pca.plotData$s, '_'))[,2]
salmon1.pca.plotData$subject  <- do.call(rbind, strsplit(salmon1.pca.plotData$s, '_'))[,3]
salmon1.pca.plotData$subject  <- gsub('don', '', salmon1.pca.plotData$subject)
salmon1.pca.plotData          <- salmon1.pca.plotData[mixedorder(as.character(salmon1.pca.plotData$day)),]
salmon1.pca.plotData$day      <- factor(salmon1.pca.plotData$day, levels = unique(salmon1.pca.plotData$day))
salmon1.pca.plotData$genotype <- factor(salmon1.pca.plotData$genotype, levels = c('none', 'KD', 'WT', 'Y664F'))

report$salmon1.pca.plot <- 
  ggplot(salmon1.pca.plotData, aes(x=x, y=y, color = genotype, shape = day, label = subject)) +
  theme_bw() +
  geom_point(size = 4, stroke = 1.5) +
  scale_shape_manual(values = 21:25) +
  geom_text(size = 5, nudge_x = 5, nudge_y = 5, show.legend = FALSE) +
  scale_color_manual(values=c('red', 'blue', 'green4', 'black')) +
  labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
       y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)'))


# Alternative plots
salmon1.rlogMat.noY664F <- salmon1.rlogMat[, - grep('Y664F', colnames(salmon1.rlogMat))]
salmon1.noY664F.pca <- prcomp(t(salmon1.rlogMat.noY664F), scale = FALSE, center = TRUE)

salmon1.pca.plotData2          <- data.frame(s = row.names(salmon1.noY664F.pca$x), x = salmon1.noY664F.pca$x[,1], y = salmon1.noY664F.pca$x[,2], z = salmon1.noY664F.pca$x[,3])
salmon1.pca.plotData2$genotype <- do.call(rbind, strsplit(salmon1.pca.plotData2$s, '_'))[,1]
salmon1.pca.plotData2$day      <- do.call(rbind, strsplit(salmon1.pca.plotData2$s, '_'))[,2]
salmon1.pca.plotData2$subject  <- do.call(rbind, strsplit(salmon1.pca.plotData2$s, '_'))[,3]
salmon1.pca.plotData2$subject  <- gsub('don', '', salmon1.pca.plotData2$subject)
salmon1.pca.plotData2          <- salmon1.pca.plotData2[mixedorder(as.character(salmon1.pca.plotData2$day)),]
salmon1.pca.plotData2$day      <- factor(salmon1.pca.plotData2$day, levels = unique(salmon1.pca.plotData2$day))
salmon1.pca.plotData2$genotype <- factor(salmon1.pca.plotData2$genotype, levels = c('none', 'KD', 'WT'))

report$salmon1.pca.plot2 <- 
  ggplot(salmon1.pca.plotData2, aes(x=x, y=y, color = genotype, shape = day, label = subject)) +
  theme_bw() +
  geom_point(size = 4, stroke = 1.5) +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values=c('black', 'red', 'dodgerblue2')) +
  labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
       y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


report$salmon1.pca.plot3 <- 
  ggplot(salmon1.pca.plotData, aes(x=x, y=y, color = genotype, shape = day, label = subject)) +
  theme_bw() +
  geom_point(size = 4, stroke = 1.5) +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values=c('black', 'red', 'dodgerblue2', 'green4')) +
  labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
       y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



# Second SALMON run PCA.

# Transform the DESeq dataset into fold change values for each transcript.
salmon2.rld <- rlog(salmon2.ddsTxi)    # Transform count data with regularized logarithm which returns log2 data normalized to library size.
salmon2.rlogMat <- assay(salmon2.rld)  # assay() extracts the matrix of normalized values.


# Perform principle component analysis.
i <- unname(apply(salmon2.rlogMat, 1, function(x) all(x < 14)))
salmon2.pca <- prcomp(t(salmon2.rlogMat[i,]), scale = FALSE, center = TRUE)

# Use the 'genotype_timePoint_donor' formatted data points to extract data for the plots
salmon2.pca.plotData          <- data.frame(s = row.names(salmon2.pca$x), x = salmon2.pca$x[,1], y = salmon2.pca$x[,2], z = salmon2.pca$x[,3])
salmon2.pca.plotData$genotype <- do.call(rbind, strsplit(salmon2.pca.plotData$s, '_'))[,1]
salmon2.pca.plotData$day      <- do.call(rbind, strsplit(salmon2.pca.plotData$s, '_'))[,2]
salmon2.pca.plotData$subject  <- do.call(rbind, strsplit(salmon2.pca.plotData$s, '_'))[,3]
salmon2.pca.plotData$subject  <- gsub('don', '', salmon2.pca.plotData$subject)
salmon2.pca.plotData          <- salmon2.pca.plotData[mixedorder(as.character(salmon2.pca.plotData$day)),]
salmon2.pca.plotData$day      <- factor(salmon2.pca.plotData$day, levels = unique(salmon2.pca.plotData$day))
salmon2.pca.plotData$genotype <- factor(salmon2.pca.plotData$genotype, levels = c('KD', 'WT', 'TrpM'))

report$salmon2.pca.plot <- 
  ggplot(salmon2.pca.plotData, aes(x=x, y=y, color = genotype, shape = day, label = subject)) +
  theme_bw() +
  geom_point(size = 4, stroke = 1.5) +
  scale_shape_manual(values = 21:25) +
  geom_text(size = 5, nudge_x = 3, nudge_y = 3, show.legend = FALSE) +
  scale_color_manual(values=c('blue', 'green4', 'darkorange')) +
  labs(x = paste0('PC1 (', sprintf("%.2f", summary(salmon1.pca)$importance[3,][1] * 100), '%)'),
       y = paste0('PC2 (', sprintf("%.2f", (summary(salmon1.pca)$importance[3,][2] - summary(salmon1.pca)$importance[3,][1])  * 100), '%)'))



# Create tx2gene lookup table from all result files.
# This will be used to add gene names to different contrasts.
d <- bind_rows(lapply(c(salmon1.files, salmon2.files), function(f){
  read.table(f, sep = '\t', header = TRUE, quote = '', comment.char = '')
}))

d$Name <- as.character(d$Name)
d <- d[!duplicated(d$Name),]
tx2gene <- data.frame(transcript_id = d$Name, 
                      gene_id       = unlist(lapply(strsplit(d$Name, '\\|'), '[[', 2)), 
                      gene_name     = unlist(lapply(strsplit(d$Name, '\\|'), '[[', 6)))


# Convenience function to add gene names from tx2gene to contrasts by using transcript id row names as a lookup.
addGeneNamesToContrast <- function(contrast){
  contrast$gene <- tx2gene[match(row.names(contrast), tx2gene$gene_id),]$gene_name
  contrast
}


# Create contrasts which contain the log2 fold change values and adjusted pvalues.

RNAseq1_Y664F_D9_vs_WT_D9       <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d9", "WT_d9")))
RNAseq1_Y664F_D12_vs_WT_D12     <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d12", "WT_d12")))
RNAseq1_Y664F_D9_vs_KD_D9       <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "Y664F_d9", "KD_d9")))
RNAseq1_WT_D9_vs_KD_D9          <- addGeneNamesToContrast(results(salmon1.ddsTxi, contrast=c("repGrp", "WT_d9", "KD_d9")))
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


# Inspect select genes to see relationships between counts and fold changes / p-values.
plotCounts(salmon2.ddsTxi, gene="ENSG00000003147.18", intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 ICA1    padj 1e-06  fold change:  2
plotCounts(salmon2.ddsTxi, gene="ENSG00000063046.17", intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 EIF4B   padj 2e-20  fold change: -2
plotCounts(salmon2.ddsTxi, gene="ENSG00000165029.15", intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 ABCA1   padj 8e-230 fold change:  7
plotCounts(salmon2.ddsTxi, gene="ENSG00000080503.23", intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 SMARCA2 padj 5e-88  fold change: -2
plotCounts(salmon2.ddsTxi, gene="ENSG00000077080.9",  intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 ACTL6B  padj NA     fold change: 0.9
plotCounts(salmon2.ddsTxi, gene="ENSG00000117713.19", intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 ARID1A  padj 0.001  fold change: -1

# Big fold change from minority of samples showing large changes.
plotCounts(salmon1.ddsTxi, gene="ENSG00000277067.4", intgroup="repGrp")  # RNAseq2_WT_D9_vs_KD_D9 ARID1A  padj 2.973e-07  fold change: 39

save.image(file='savePoints/sp1.RData')



# Heatmaps
#--------------------------------------------------------------------------------------------------

RNAseq1_Y664F_vs_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_Y664F_D9_vs_KD_D9),   exp = 'Y664F D9'),
                   dplyr::mutate(data.frame(RNAseq1_Y664F_D12_vs_KD_D12), exp = 'Y664F D12')) %>%
  dplyr::mutate(exp = factor(exp, levels = c('Y664F D9', 'Y664F D12'))) %>%
  dplyr::filter(! is.na(padj)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup()



# RNAseq2  WT / TrpM vs KD

RNAseq2_WT_TrpM_vs_KD <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq2_WT_D6_vs_KD_D6),   exp = 'WT D6'),
                   dplyr::mutate(data.frame(RNAseq2_TrpM_D6_vs_KD_D6), exp = 'TrpM D6'),
                   dplyr::mutate(data.frame(RNAseq2_WT_D9_vs_KD_D9),   exp = 'WT D9'),
                   dplyr::mutate(data.frame(RNAseq2_TrpM_D9_vs_KD_D9), exp = 'TrpM D9')) %>%
  dplyr::mutate(exp = factor(exp, levels = c('WT D6', 'TrpM D6', 'WT D9', 'TrpM D9'))) %>%
  dplyr::filter(! is.na(padj)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup()

RNAseq2_WT_TrpM_vs_KD.genes <- 
  dplyr::group_by(RNAseq2_WT_TrpM_vs_KD, gene) %>%
  dplyr::mutate(exps = n_distinct(exp)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(! is.na(padj) & abs(log2FoldChange) >= 3 & abs(log2FoldChange) <= 14 & exps == n_distinct(exp) & padj <= 1e-70) %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  unique()

RNAseq2_WT_TrpM_vs_KD.genesForPathways <- 
  dplyr::group_by(RNAseq2_WT_TrpM_vs_KD, gene) %>%
  dplyr::mutate(exps = n_distinct(exp)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(! is.na(padj) & abs(log2FoldChange) >= 3 & abs(log2FoldChange) <= 14 & exps == n_distinct(exp) & padj <= 1e-20) %>%
  dplyr::arrange(padj)


plot.data.A <- dplyr::filter(RNAseq2_WT_TrpM_vs_KD, gene %in% RNAseq2_WT_TrpM_vs_KD.genes[1:38])  %>% dplyr::arrange(padj) %>% dplyr::mutate(gene = factor(gene, levels = rev(unique(gene)))) 
plot.data.B <- dplyr::filter(RNAseq2_WT_TrpM_vs_KD, gene %in% RNAseq2_WT_TrpM_vs_KD.genes[39:75]) %>% dplyr::arrange(padj) %>% dplyr::mutate(gene = factor(gene, levels = rev(unique(gene)))) 

report$RNAseq2_WT_TrpM_vs_KD.plotA <- 
  make_square(ggplot(plot.data.A, aes(x = exp, y = gene, fill = log2FoldChange)) +
              theme_bw() +
              geom_tile(color="gray50",size=0.6) + 
              scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", midpoint=0, limits=c(-12, 12), breaks = c(-12, -6, 0, 6, 12)) + 
              scale_y_discrete(expand=c(0,0)) +
              scale_x_discrete(expand=c(0,0)) +
              labs(x='', y = '') +
              theme(axis.text.x = element_text(angle = 90, hjust = 1),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"),
                    legend.position="bottom") +
              guides(fill=guide_colorbar(title.position = "top", barwidth=5)), fudge = 0.5)

report$RNAseq2_WT_TrpM_vs_KD.plotB <- 
  make_square(ggplot(plot.data.B, aes(x = exp, y = gene, fill = log2FoldChange)) +
                theme_bw() +
                geom_tile(color="gray50",size=0.6) + 
                scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", midpoint=0, limits=c(-12, 12), breaks = c(-12, -6, 0, 6, 12)) + 
                scale_y_discrete(expand=c(0,0)) +
                scale_x_discrete(expand=c(0,0)) +
                labs(x='', y = '') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black"),
                      legend.position="bottom") +
                guides(fill=guide_colorbar(title.position = "top", barwidth=5)), fudge = 0.5)



# Heat maps for select gene lists.
#--------------------------------------------------------------------------------------------------

f <- list()
f[['Signaling molecules']] <- 'data/Pawlicki_signaling.txt'
f[['TCR a/b']] <- 'data/Pawlicki_TCRab.txt'
f[['Transcriptional regulators']] <- 'data/Pawlicki_TranscriptionalRegulators.txt'
f[['Upregulated in WNT']] <- 'data/Pawlicki_upInWT.txt'

g <- bind_rows(mapply(function(n, x){
       genes <- readLines(x)
       o <- bind_rows(subset(RNAseq2_WT_TrpM_vs_KD, gene %in% genes),
                      subset(RNAseq1_Y664F_vs_KD, gene %in% genes))
       
       o$title <- n
       o
}, names(f), f, SIMPLIFY = FALSE))

g$exp <- factor(as.character(g$exp), levels = c("WT D6", "TrpM D6", "WT D9", "TrpM D9", "Y664F D9", "Y664F D12"))
g$marker <- ifelse(g$padj <= 0.05, TRUE, FALSE)

options(digits=1)
report$genesOfInterest_RNAseq <-
  ggplot(g, aes(x = exp, y = gene, fill = log2FoldChange, shape = marker)) +
         theme_bw() +
         geom_tile(color = 'gray20') +
         geom_point(alpha = 0.25) +
         scale_shape_manual(values = c(1, 32)) +
         scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", 
                              midpoint=0, limits=c(-11.1, 11.1), breaks = c(-11, -5, 0, 5, 11)) +
         scale_y_discrete(expand=c(0,0)) +
         scale_x_discrete(expand=c(0,0)) +
         labs(x='', y = '') +
         theme(strip.background = element_blank(),
               axis.text.x = element_text(angle = 90, hjust = 1),
               panel.border = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               legend.position="bottom") +
        guides(shape = FALSE, fill=guide_colorbar(title.position = "top",  barwidth=10)) +
        facet_wrap(~title, scales = 'free') 



# SWI/SNF
#--------------------------------------------------------------------------------------------------
genes <- RNAseq2_WT_D9_vs_KD_D9$gene[grepl('^SMARC|^ARID1|^ACTL6', RNAseq2_WT_D9_vs_KD_D9$gene)]

g <- bind_rows(subset(RNAseq2_WT_TrpM_vs_KD, gene %in% genes),
               subset(RNAseq1_Y664F_vs_KD, gene %in% genes))
g$exp <- factor(as.character(g$exp), levels = c("WT D6",  "WT D9", "TrpM D6", "TrpM D9", "Y664F D9", "Y664F D12"))

g$exp <- factor(as.character(g$exp), levels = c("WT D6", "TrpM D6", "WT D9", "TrpM D9", "Y664F D9", "Y664F D12"))
g$marker <- ifelse(g$padj <= 0.05, TRUE, FALSE)

options(digits=1)
report$SWI_SNF_genes_RNAseq <-
  ggplot(g, aes(x = exp, y = gene, fill = log2FoldChange, shape = marker)) +
  theme_bw() +
  geom_tile(color = 'gray20') +
  geom_point(alpha = 0.25) +
  scale_shape_manual(values = c(1, 32)) +
  scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(-6.1, 6.1), breaks = c(-6, -3, 0, 3, 6)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x='', y = '') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="bottom") +
  guides(shape = FALSE, fill=guide_colorbar(title.position = "top",  barwidth=10)) 





# Volcano plots
#--------------------------------------------------------------------------------------------------

createVolcanoPlot <- function(RNAseq, title, file){
  plot.data <- subset(data.frame(RNAseq), ! is.na(padj))
  plot.data$class <- ifelse(plot.data$padj > 0.05, 'No significant change', ifelse(plot.data$log2FoldChange >= 0, 'Increased', 'Decreased'))
  plot.data$class <- factor(as.character(plot.data$class), levels = c('No significant change', 'Increased', 'Decreased'))

  plot.data$foldChange <- 2^plot.data$log2FoldChange
  gc()
  write.xlsx2(plot.data, file = sub('\\.\\S+$', '.xlsx', file), col.names = TRUE, row.names = TRUE)
  
  p <- ggplot(plot.data, aes(foldChange, -log10(padj), fill = class)) +
       theme_bw() +
       scale_fill_manual(name = 'Transcription', values = c('gray50', 'red', 'dodgerblue2')) +
       geom_point(shape = 21, size = 2, alpha = 0.5, color = 'black') +
       scale_x_log10() +
       annotation_logticks(base = 10, sides="b") +
       theme(panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             axis.line = element_line(colour = "black")) +
      ggtitle(title) +
      guides(fill = guide_legend(override.aes=list(shape=21, size = 3)),
             shape = guide_legend(override.aes=list(size = 3))) +
      labs(x= 'Fold Change', y = '-log10(adjusted p-value)')
  
  ggsave(p, file = file)
  p
}

report$RNAseq1_WT_D6_vs_KD_D6_volcano   <- createVolcanoPlot(RNAseq1_WT_D6_vs_KD_D6,   'WT D6 vs KD D6',  'paper_figures_and_tables/RNAseq1_WT_D6_vs_KD_D6_volcano.pdf')
report$RNAseq1_WT_D9_vs_KD_D9_volcano   <- createVolcanoPlot(RNAseq1_WT_D9_vs_KD_D9,   'WT D9 vs KD D9',  'paper_figures_and_tables/RNAseq1_WT_D9_vs_KD_D9_volcano.pdf')
report$RNAseq1_WT_D12_vs_KD_D12_volcano <- createVolcanoPlot(RNAseq1_WT_D12_vs_KD_D12, 'WT D12 vs KD D12', 'paper_figures_and_tables/RNAseq1_WT_D12_vs_KD_D12_volcano.pdf')

report$RNAseq1_WT_D33_vs_KD_D12_volcano <- createVolcanoPlot(RNAseq1_WT_D33_vs_KD_D12, 'WT D33 vs KD D12', 'paper_figures_and_tables/RNAseq1_WT_D33_vs_KD_D12_volcano.pdf')
report$RNAseq1_WT_D63_vs_KD_D12_volcano <- createVolcanoPlot(RNAseq1_WT_D63_vs_KD_D12, 'WT D63 vs KD D12', 'paper_figures_and_tables/RNAseq1_WT_D63_vs_KD_D12_volcano.pdf')


#--------------------------------------------------------------------------------------------------

RNAseq2_TrpM_vs_WT <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq2_TrpM_D6_vs_WT_D6), exp = 'D6'),
                   dplyr::mutate(data.frame(RNAseq2_TrpM_D9_vs_WT_D9), exp = 'D9')) %>%
     dplyr::mutate(exp = factor(exp, levels = c('D6', 'D9'))) %>%
     dplyr::filter(! is.na(padj)) %>%
     dplyr::group_by(exp, gene) %>% 
     dplyr::top_n(-1, wt = padj) %>%
     dplyr::ungroup()

RNAseq2_TrpM_vs_WT.genes <- 
  dplyr::group_by(RNAseq2_TrpM_vs_WT, gene) %>%
  dplyr::mutate(exps = n_distinct(exp)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(! is.na(padj) & abs(log2FoldChange) >= 4 & abs(log2FoldChange) <= 14 & exps == n_distinct(exp) & padj <= 1e-15) %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  unique()

plot.data <- dplyr::filter(RNAseq2_TrpM_vs_WT, gene %in% RNAseq2_TrpM_vs_WT.genes) %>% dplyr::arrange(padj) %>% dplyr::mutate(gene = factor(gene, levels = rev(unique(gene)))) 


report$RNAseq2.pval.TrpM_vs_WT <- 
  make_square(ggplot(plot.data, aes(x = exp, y = gene, fill = log2FoldChange)) +
                theme_bw() +
                geom_tile(color="gray50",size=0.6) + 
                scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", midpoint=0, limits=c(-8, 8), breaks = c(-8, -4, 0, 4, 8)) + 
                scale_y_discrete(expand=c(0,0)) +
                scale_x_discrete(expand=c(0,0)) +
                labs(x='', y = '') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black"),
                      legend.position="bottom") +
                guides(fill=guide_colorbar(title.position = "top", barwidth=5)), fudge = 0.5)


#----------------------




#--------------------------------------------------------------------------------------------------
# RNAseq1_Y664F_D9_vs_WT_D9

RNAseq1_Y664F_vs_WT <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_Y664F_D9_vs_WT_D9),   exp = 'Y664F D9'),
                   dplyr::mutate(data.frame(RNAseq1_Y664F_D12_vs_WT_D12), exp = 'Y664F D12')) %>%
     dplyr::mutate(exp = factor(exp, levels = c('Y664F D9', 'Y664F D12'))) %>%
     dplyr::filter(! is.na(padj)) %>%
     dplyr::group_by(exp, gene) %>% 
     dplyr::top_n(-1, wt = padj) %>%
     dplyr::ungroup()

RNAseq1_Y664F_vs_WT.genes <- 
  dplyr::group_by(RNAseq1_Y664F_vs_WT, gene) %>%
  dplyr::mutate(exps = n_distinct(as.character(exp))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(abs(log2FoldChange) >= 2 & abs(log2FoldChange) <= 14 & exps == n_distinct(exp) & padj <= 1e-03) %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  unique()

plot.data <- dplyr::filter(RNAseq1_Y664F_vs_WT, gene %in% RNAseq1_Y664F_vs_WT.genes) %>% dplyr::arrange(padj) %>% dplyr::mutate(gene = factor(gene, levels = rev(unique(gene)))) 

report$RNAseq1.pval.WT_Y664F_vs_WT <- 
  make_square(ggplot(plot.data, aes(x = exp, y = gene, fill = log2FoldChange)) +
                theme_bw() +
                geom_tile(color="gray50",size=0.6) + 
                scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", midpoint=0, limits=c(-10, 10), breaks = c(-10, -5, 0, 5, 10)) + 
                scale_y_discrete(expand=c(0,0)) +
                scale_x_discrete(expand=c(0,0)) +
                labs(x='', y = '') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black"),
                      legend.position="bottom") +
                guides(fill=guide_colorbar(title.position = "top", barwidth=5)), fudge = 0.5)


#--------------------------------------------------------------------------------------------------


RNAseq1_WT_vs_earlyWT <- 
  dplyr::bind_rows(dplyr::mutate(data.frame(RNAseq1_WT_D9_vs_WT_D6),  exp = 'D9'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D12_vs_WT_D6), exp = 'D12'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D33_vs_WT_D6), exp = 'D33'),
                   dplyr::mutate(data.frame(RNAseq1_WT_D63_vs_WT_D6), exp = 'D63')) %>%
  dplyr::mutate(exp = factor(exp, levels = c('D9', 'D12', 'D33', 'D63'))) %>%
  dplyr::filter(! is.na(padj)) %>%
  dplyr::group_by(exp, gene) %>% 
  dplyr::top_n(-1, wt = padj) %>%
  dplyr::ungroup()

RNAseq1_WT_vs_earlyWT.genes <- subset(data.frame(RNAseq1_WT_D63_vs_WT_D6), padj <= 1e-25 & abs(log2FoldChange) >= 3 & abs(log2FoldChange) <= 14)$gene

RNAseq1_WT_vs_earlyWT.genesForPathways <- subset(data.frame(RNAseq1_WT_D63_vs_WT_D6), padj <= 1e-10 & abs(log2FoldChange) >= 3 & abs(log2FoldChange) <= 14)$gene


plot.data <- dplyr::filter(RNAseq1_WT_vs_earlyWT, gene %in% RNAseq1_WT_vs_earlyWT.genes) %>% dplyr::arrange(padj) %>% dplyr::mutate(gene = factor(gene, levels = rev(unique(gene)))) 

report$RNAseq1.pval.WT_vs_earlyWT <- 
  make_square(ggplot(plot.data, aes(x = exp, y = gene, fill = log2FoldChange)) +
                theme_bw() +
                geom_tile(color="gray50",size=0.6) + 
                scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", midpoint=0, limits=c(-8, 8), breaks = c(-8, -4, 0, 8, 4)) + 
                scale_y_discrete(expand=c(0,0)) +
                scale_x_discrete(expand=c(0,0)) +
                labs(x='', y = '') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black"),
                      legend.position="bottom") +
                guides(fill=guide_colorbar(title.position = "top", barwidth=5)), fudge = 0.5)


save.image(file='savePoints/sp2.RData')





#--------------------------------------------------------------------------------------------------



# Add entrez ids to the RNAseq data.

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

### mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), mart = ensembl)
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "entrezgene_accession"), mart = ensembl)
mapping$entrezgene <- mapping$entrezgene_id

RNAseq1_WT_vs_earlyWT$ensembl     <- sub('\\.\\d+$', '', tx2gene[match(toupper(RNAseq1_WT_vs_earlyWT$gene),     toupper(tx2gene$gene_name)),]$gene_id)
RNAseq1_Y664F_D9_vs_WT_D9$ensembl <- sub('\\.\\d+$', '', tx2gene[match(toupper(RNAseq1_Y664F_D9_vs_WT_D9$gene), toupper(tx2gene$gene_name)),]$gene_id)
RNAseq2_WT_TrpM_vs_KD$ensembl     <- sub('\\.\\d+$', '', tx2gene[match(toupper(RNAseq2_WT_TrpM_vs_KD$gene),     toupper(tx2gene$gene_name)),]$gene_id)


RNAseq1_WT_vs_earlyWT     <- subset(RNAseq1_WT_vs_earlyWT,     abs(log2FoldChange) >= 2 & padj <= 1e-3)
RNAseq1_Y664F_D9_vs_WT_D9 <- subset(RNAseq1_Y664F_D9_vs_WT_D9, abs(log2FoldChange) >= 2 & padj <= 1e-3)
RNAseq2_WT_TrpM_vs_KD     <- subset(RNAseq2_WT_TrpM_vs_KD,     abs(log2FoldChange) >= 2 & padj <= 1e-3)


string_db <- STRINGdb$new( version = "10", species = 9606, score_threshold = 0, input_directory = report$STRINGdb_dataFiles)

RNAseq1_WT_vs_earlyWT <- string_db$map(data.frame(RNAseq1_WT_vs_earlyWT), "gene", removeUnmappedRows = FALSE)
i <- which(is.na(RNAseq1_WT_vs_earlyWT$STRING_id))
RNAseq1_WT_vs_earlyWT[i,]$STRING_id <- report$string_db.alt.aliases[match(toupper(RNAseq1_WT_vs_earlyWT[i,]$gene), toupper(report$string_db.alt.aliases$ids)),]$STRING_id

RNAseq1_Y664F_D9_vs_WT_D9 <- string_db$map(data.frame(RNAseq1_Y664F_D9_vs_WT_D9), "gene", removeUnmappedRows = FALSE)
i <- which(is.na(RNAseq1_Y664F_D9_vs_WT_D9$STRING_id))
RNAseq1_Y664F_D9_vs_WT_D9[i,]$STRING_id <- report$string_db.alt.aliases[match(toupper(RNAseq1_Y664F_D9_vs_WT_D9[i,]$gene), toupper(report$string_db.alt.aliases$ids)),]$STRING_id

RNAseq2_WT_TrpM_vs_KD <- string_db$map(data.frame(RNAseq2_WT_TrpM_vs_KD), "gene", removeUnmappedRows = FALSE)
i <- which(is.na(RNAseq2_WT_TrpM_vs_KD$STRING_id))
RNAseq2_WT_TrpM_vs_KD[i,]$STRING_id <- report$string_db.alt.aliases[match(toupper(RNAseq2_WT_TrpM_vs_KD[i,]$gene), toupper(report$string_db.alt.aliases$ids)),]$STRING_id


RNAseq1_WT_vs_earlyWT$entrezgene     <- mapping[match(RNAseq1_WT_vs_earlyWT$ensembl, mapping$ensembl_gene_id),]$entrezgene
RNAseq1_Y664F_D9_vs_WT_D9$entrezgene <- mapping[match(RNAseq1_Y664F_D9_vs_WT_D9$ensembl, mapping$ensembl_gene_id),]$entrezgene
RNAseq2_WT_TrpM_vs_KD$entrezgene     <- mapping[match(RNAseq2_WT_TrpM_vs_KD$ensembl, mapping$ensembl_gene_id),]$entrezgene


# KEGG pathway enrichments
#--------------------------------------------------------------------------------------------------
# Detach the RMySQL package because it interferes with STRINGdb
detach('package:RMySQL', unload = TRUE, character.only = TRUE)

# Function which creates a table of enriched KEGG terms from associated STRINGdb protein ids.
# The names of genes considered in the enrichments are provided and are listed with uppercase 
# letters if the gene showed increased transcription levels and listed with lowercase letters 
# if the gene showed decreased transcription levels.

createKEGGenrichmentTable <- function(o){
  e <- string_db$get_enrichment(o$STRING_id, category = 'KEGG', methodMT = "fdr", iea = FALSE)
  e <- subset(e, pvalue_fdr <= 0.05)
  e$geneListLength <- n_distinct(o$STRING_id)
  
  p <- string_db$get_term_proteins(e$term_id)
  p <- p[p$STRING_id %in% o$STRING_id,]
  
  o.up   <- subset(o, log2FoldChange >= 0)$STRING_id
  o.down <- subset(o, log2FoldChange < 0)$STRING_id
  
  p <- dplyr::rowwise(p) %>% dplyr::mutate(preferred_name2 = ifelse(STRING_id %in% o.up, toupper(preferred_name), tolower(preferred_name))) %>% dplyr::ungroup()
  p <- dplyr::group_by(p, term_id) %>% dplyr::summarise(genes = paste0(sort(preferred_name2), collapse = ', ')) %>% dplyr::ungroup()
  dplyr::left_join(e, p, by = 'term_id')
}


# Function which creates KEGG schematics color coded by gene trascription fold changes.

createKEGGschematics <- function(o, schematicPrefix = 'exp', pathways = NA){
  foldchanges = unname(o$log2FoldChange)
  names(foldchanges) = o$entrezgene
  #browser()
  pathview(gene.data=foldchanges, pathway.id = pathways, species="hsa", out.suffix = schematicPrefix)

  if(! dir.exists('KEGG.schematics')) dir.create('KEGG.schematics')
  system(paste('mv', paste(list.files(pattern = schematicPrefix), collapse = ' '), 'KEGG.schematics/'))
  unlink(list.files(pattern = paste0('hsa', pathways, collapse = '|')))
}


# WT_D9_vs_KD_D9 -- KEGG enrichment analysis with STRINGdb.  
report$WT_D9_vs_KD_D9.KEGG <- createKEGGenrichmentTable(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'WT D9' & 
                                                                 ! is.na(log2FoldChange) &  
                                                                 ! is.na(STRING_id) & 
                                                                 padj <= 1e-5 & 
                                                                 abs(log2FoldChange) >= 2 &
                                                                 abs(log2FoldChange) <= 14))

createKEGGschematics(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'WT D9' & 
                              ! is.na(log2FoldChange) &  
                              ! is.na(entrezgene) & 
                              padj <= 1e-5 & 
                              abs(log2FoldChange) >= 2 &
                              abs(log2FoldChange) <= 14),
                     schematicPrefix = 'WT_D9_vs_KD_D9',
                     pathways = report$WT_D9_vs_KD_D9.KEGG[1:10,]$term_id)



# 
report$TrpM_D9_vs_KD_D9.KEGG <- createKEGGenrichmentTable(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'TrpM D9' & 
                                                                   ! is.na(log2FoldChange) &  
                                                                   ! is.na(STRING_id) & 
                                                                   padj <= 1e-5 & 
                                                                   abs(log2FoldChange) >= 2 &
                                                                   abs(log2FoldChange) <= 14))


createKEGGschematics(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'TrpM D9' & 
                              ! is.na(log2FoldChange) &  
                              ! is.na(STRING_id) & 
                              padj <= 1e-5 & 
                              abs(log2FoldChange) >= 2 &
                              abs(log2FoldChange) <= 14),
                     schematicPrefix = 'TrpM_D9_vs_KD_D9',
                     pathways = report$TrpM_D9_vs_KD_D9.KEGG[1:10,]$term_id)




report$WT_late_vs_WT_early.KEGG <- createKEGGenrichmentTable(subset(RNAseq1_WT_vs_earlyWT, exp == 'D63' & 
                                                                   ! is.na(log2FoldChange) &  
                                                                   ! is.na(STRING_id) & 
                                                                   padj <= 1e-5 & 
                                                                   abs(log2FoldChange) >= 2 &
                                                                   abs(log2FoldChange) <= 14))

createKEGGschematics(subset(RNAseq1_WT_vs_earlyWT, exp == 'D63' & 
                              ! is.na(log2FoldChange) &  
                              ! is.na(STRING_id) & 
                              padj <= 1e-5 & 
                              abs(log2FoldChange) >= 2 &
                              abs(log2FoldChange) <= 14),
                     schematicPrefix = 'WT_late_vs_WT_early',
                     pathways = report$WT_late_vs_WT_early.KEGG[1:10,]$term_id)



# GO term enrichment
#--------------------------------------------------------------------------------------------------
report$WT_D9_vs_KD_D9.GO <- 
  string_db$get_enrichment(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'WT D9' & 
                                    ! is.na(log2FoldChange) &  
                                    ! is.na(STRING_id) & 
                                    padj <= 1e-5 & 
                                    abs(log2FoldChange) >= 2 &
                                    abs(log2FoldChange) <= 14)$STRING_id, category = 'Process', methodMT = "fdr", iea = FALSE)
report$WT_D9_vs_KD_D9.GO <- subset(report$WT_D9_vs_KD_D9.GO, pvalue_fdr <= 0.05)


report$TrpM_D9_vs_KD_D9.GO <- 
  string_db$get_enrichment(subset(RNAseq2_WT_TrpM_vs_KD, exp == 'TrpM D9' & 
                                    ! is.na(log2FoldChange) &  
                                    ! is.na(STRING_id) & 
                                    padj <= 1e-5 & 
                                    abs(log2FoldChange) >= 2 &
                                    abs(log2FoldChange) <= 14)$STRING_id, category = 'Process', methodMT = "fdr", iea = FALSE)
report$TrpM_D9_vs_KD_D9.GO <- subset(report$TrpM_D9_vs_KD_D9.GO , pvalue_fdr <= 0.05)



save.image(file='savePoints/sp3.RData')



# Integration site data
#-------------------------------------------------------------------------------------------------------

# Read in sample data.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="NPM_ALK"')
dbDisconnect(dbConn)

# Retrieve and standardize fragments, call intSites, calculate abundances, and annotate sites.
report$intSites <- getDBgenomicFragments(samples = samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
                   stdIntSiteFragments() %>%
                   collapseReplicatesCalcAbunds() %>%
                   annotateIntSites()


save.image(file = 'savePoints/sp3.RData')


report$intSites$genotype <- ifelse(grepl('Y664F', report$intSites$patient), 'Y664F', 'WT')


# Create a list of relative abundance plots for each subject.
library(RColorBrewer)
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
    theme(plot.title = element_text(size = 8)) +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    scale_y_continuous(labels = scales::percent)
})


# Detach the RMySQL package because it interferes with STRINGdb
detach('package:RMySQL', unload = TRUE, character.only = TRUE)


# Create list of patients with multiple time points.
report$subjectsWithMultTimePoints <- 
  dplyr::group_by(data.frame(report$intSites), patient) %>% 
  dplyr::summarise(tps = n_distinct(timePoint)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(tps > 1) %>% 
  dplyr::select(patient) %>% 
  unlist() %>% unname()



# Map STRINGdb ids 
#string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory = report$STRINGdb_dataFiles)

d <- string_db$map(data.frame(report$intSites), "nearestFeature", removeUnmappedRows = FALSE)
i <- which(is.na(d$STRING_id));  message(paste0('gene ids not mapped to STRINGdb: ', round((length(i) / nrow(d))*100, 2), '%'))
d[i,]$STRING_id <- report$string_db.alt.aliases[match(toupper(d[i,]$nearestFeature), toupper(report$string_db.alt.aliases$ids)),]$STRING_id
i <- which(is.na(d$STRING_id));  message(paste0('gene ids not mapped to STRINGdb: ', round((length(i) / nrow(d))*100, 2), '%'))
report$intSites$STRING_id <- d$STRING_id



# Create an integration frequency plot for WT subjects.
library(ggrepel)
report$preVspostFreqplot.WT <- preVsPostFreqPlot(subset(report$intSites, patient %in% c('pNA92', 'pNA93', 'pNA85')),
                                                 14, 
                                                 readRDS('data/humanOncoGenesList.rds'), 
                                                 nGenesToLabel = 5,
                                                 mustLabel = c('ARID1A'))


# Create an integration frequency plot for Y664F subjects.
report$preVspostFreqplot.Y664F <- preVsPostFreqPlot(subset(report$intSites, patient %in% c('pNA92_Y664F', 'pNA93_Y664F', 'pNA85_Y664F')),
                                                    14, 
                                                    readRDS('data/humanOncoGenesList.rds'), 
                                                    nGenesToLabel = 5,
                                                    mustLabel = c('ARID1A'))



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



#--------------------------------------------------------------------------------------------------


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

save.image(file = 'savePoints/sp4.RData')
saveRDS(report, file = 'report.rds')

library(GenomicRanges)
library(STRINGdb)
library(biomaRt)
library(DESeq2)
library(gtools)
library(tidyverse)
library(pathview)
library(ggrepel)
library(RColorBrewer)
library(fgsea)
source('./lib.R')

load('GSEA.img')


create_GSEA_hallmark_plots <- function(d, hallmark, label){
  invisible(lapply(unique(d$timePoint), function(tp){
    o <- subset(d, timePoint == tp)
    ranks <- o$stat
    names(ranks) <- o$entrezgene
    ranks <- ranks[! is.na(names(ranks))]
    
    p <- ggplot(data.frame(y = sort(ranks, decreasing = TRUE), x = 1:length(ranks)), aes(x, y)) + 
      theme_bw() +
      geom_bar(stat = "identity", width=1.05, fill="black") +
      scale_x_continuous(label=scales::comma) +
      theme(axis.text=element_text(size=12),  axis.title = element_text(size=14), plot.title = element_text(size = 16), plot.margin = unit(c(1,1,1,1), "cm")) +
      ggtitle(paste0(tp, ' ordered gene list stats')) +
      labs(x = 'Gene', y = 'Stat')
    
    ggsave(p, file = paste0('figures_and_tables/GSEA/', label, '_', tp, '_ordered_gene_list_histogram.pdf'))
    
    invisible(mapply(function(pathway, proteins){
      set.seed(46)
      x <- plotEnrichment(proteins, ranks, gseaParam = 1)
      pathwayList <- list()
      pathwayList[[pathway]] <- proteins
      fgseaRes <- fgsea(pathwayList, ranks, nperm=10000)
      p <- fgseaEnrichmentPlotFormat(x, pathway)
      ggsave(p, file = paste0('figures_and_tables/GSEA/', label, '_', tp, '_', paste0(sprintf("%.2f", fgseaRes$NES), '_', pathway, '.pdf')))
    }, names(hallmark), hallmark))
  }))
}

# WT_vs_KD found in both first and second RNAseq experiments.

create_GSEA_hallmark_plots(RNAseq1_Y664F_vs_KD, hallmark, 'RNAseq1.Y664F_vs_KD')
create_GSEA_hallmark_plots(RNAseq1_Y664F_vs_WT, hallmark, 'RNAseq1.Y664F_vs_WT')
create_GSEA_hallmark_plots(RNAseq1_WT_vs_KD,    hallmark, 'RNAseq1.WT_vs_KD')
create_GSEA_hallmark_plots(RNAseq2_WT_vs_KD,    hallmark, 'RNAseq2.WT_vs_KD')
create_GSEA_hallmark_plots(RNAseq2_TrpM_vs_WT,  hallmark, 'RNAseq2.TrpM_vs_WT')
create_GSEA_hallmark_plots(RNAseq2_TrpM_vs_KD,  hallmark, 'RNAseq2.TrpM_vs_KD')

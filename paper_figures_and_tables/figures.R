library(knitr)
library(kableExtra)
library(png)
library(grid)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(grDevices)
library(gtools)
ppNum  <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
report <- readRDS('report.rds')


f <- 
  data.frame(subset(report$intSites, report$intSites$patient %in% c('pNA92', 'pNA93', 'pNA85'))) %>%
  dplyr::group_by(patient, timePointDays) %>%
  dplyr::summarise(clones = n_distinct(posid)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(genome = ifelse(grepl('Y664F', patient), 'Y664F', 'WT')) %>%
  ggplot(aes(timePointDays, log10(clones), group = patient, color = patient)) + 
  theme_bw() + 
  scale_color_manual(name = 'Subject', values = c('red', 'blue', 'green3')) +
  geom_point(size=3) + 
  geom_line() + 
  labs(x = 'Days', y = 'log10 (Unique clones)') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        
ggsave(f, file = "WT_clonality_decrease.pdf")


x <- report$sampleAbundancePlotsData$pNA93
x <- x[rev(order(x$relAbund)),]
x$geneLabel <- paste0(x$nearestFeature, '\n', x$posid)
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
  
df$timePoint <- factor(df$timePoint, levels=mixedsort(unique(df$timePoint)))
df <- arrange(df, desc(relAbund))
df$sites     <- factor(df$sites, levels=unique(df$sites))
  
f <- ggplot() + 
theme_bw() +
geom_bar(data=df, aes(timePoint, relAbund/100, fill=sites), stat = "identity") + 
guides(fill=FALSE) +
labs(x='', y='') +
ggtitle(x$patient[1]) +
scale_fill_manual(values=c('#DCDCDC', colorRampPalette(brewer.pal(12, "Paired"))(report$sampleAbundancePlotsData_n))) +
scale_y_continuous(labels = scales::percent) +
theme(axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave(f, file = "WT_clonality_decrease_NA93_example.pdf")


ggsave(report$ALK_ALCL_heatmap_a,    file = 'ALK_ALCL_heatmap_a.pdf')
ggsave(report$ALK_ALCL_heatmap_b,    file = 'ALK_ALCL_heatmap_b.pdf')
ggsave(report$Tcell_TFs_heatmap,     file = 'Tcell_TFs_heatmap.pdf')
ggsave(report$Tcell_sig_heatmap,     file = 'Tcell_sig_heatmap.pdf')
ggsave(report$Tcell_recp_heatmap,    file = 'Tcell_recp_heatmap.pdf')
ggsave(report$NonTcell_spec_heatmap, file = 'NonTcell_spec_heatmap.pdf')
ggsave(report$EmbStemSpec_heatmap,   file = 'EmbStemSpec_heatmap.pdf')



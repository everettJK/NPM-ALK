heatmap_dims <- function(p) {
  .x <- as.character(p$mapping$x)
  .y <- as.character(p$mapping$y)
  
  .x <- .x[-grep('~', .x)]
  .y <- .y[-grep('~', .y)]
  
  ncols <- length(unique(p$data[[.x]]))
  nrows <- length(unique(p$data[[.y]]))
  return(list(ncols=ncols, nrows=nrows))
}

make_square <- function(p, fudge=1) {
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows/dims$ncols)*fudge)
}

ppNum <- function (n)  format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


# Convenience function to add gene names from tx2gene to contrasts by using transcript id row names as a lookup.
# Requires global tx2gene object.
addGeneNamesToContrast <- function(contrast){
  contrast$gene <- tx2gene[match(row.names(contrast), tx2gene$gene_id),]$gene_name
  contrast
}


createVolcanoPlot <- function(RNAseq, title, file){
  plot.data <- subset(data.frame(RNAseq), ! is.na(padj) & abs(log2FoldChange) <= 15)
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



# Convenience function for adding truncated ensembl ids to the RNAseq data, mapping STRING ids, 
# expanding STRING annotations with aliases from HGNC, and then adding entrez gene names via bioMart.
# This function is dependent upon global objects tx2gene, string_db, report$string_db.alt.aliases, and mapping.
addAndExtendSTRINGids <- function(d){
  d$ensembl       <- sub('\\.\\d+$', '', tx2gene[match(toupper(d$gene), toupper(tx2gene$gene_name)),]$gene_id)
  d <- subset(d, abs(log2FoldChange) >= 2 & padj <= 1e-3)
  
  d <- string_db$map(data.frame(d), "gene", removeUnmappedRows = FALSE)
  i <- which(is.na(d$STRING_id))
  d[i,]$STRING_id <- report$string_db.alt.aliases[match(toupper(d[i,]$gene), toupper(report$string_db.alt.aliases$ids)),]$STRING_id
  
  d$entrezgene     <- mapping[match(d$ensembl, mapping$ensembl_gene_id),]$entrezgene
  d
}





createGeneListHeatMap <- function(data, geneList, scaleLimit, file, orderByFoldChange = TRUE){
  plot.data      <- subset(data, gene %in% geneList)
  plot.data$gene <- factor(as.character(plot.data$gene), levels = rev(geneList))
  plot.data$sig  <- factor(ifelse(plot.data$padj <= 0.05, TRUE, FALSE), levels = c(TRUE, FALSE))
  plot.data$timePoint <- factor(as.character(plot.data$timePoint), levels =  gtools::mixedsort(unique(plot.data$timePoint)))
  
  midScaleMark <- floor(scaleLimit/2)
  
  if(orderByFoldChange){
    plot.data <- group_by(plot.data, gene) %>% mutate(rowsum = sum(log2FoldChange)) %>% ungroup() %>% arrange(rowsum) 
    plot.data$gene <- factor(as.character(plot.data$gene), levels = unique(plot.data$gene))
  }
  
  p <- make_square(ggplot(plot.data, aes(x = timePoint, y = gene, fill = log2FoldChange, shape = sig)) +
                theme_bw() +
                geom_tile(color="gray90",size=0.6) + 
                #geom_tile(color="black",size=0.6) +
                scale_shape_manual(values = c(42, 32)) +
                geom_point(size = 2, color = 'gray25', show.legend = FALSE) +
                scale_fill_gradient2(name = 'Fold change', low="navy", mid="white", high="red", midpoint=0, 
                                     limits=c(-scaleLimit, scaleLimit), breaks = c(-scaleLimit, -midScaleMark, 0, midScaleMark, scaleLimit)) + 
                scale_y_discrete(expand=c(0,0)) +
                scale_x_discrete(expand=c(0,0)) +
                labs(x='', y = '') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.line = element_line(colour = "black"),
                      legend.position="bottom") +
                guides(fill=guide_colorbar(title.position = "top", barwidth=5)), fudge = 0.75) 

  ggsave(p, file = file)
  p
}



timePoint2numeric <- function(x, interval = 'days'){
  timePointType <- stringr::str_match(x, '[DMY]')
  if(as.logical(is.na(timePointType))) timePointType <- 'X'
  n <- as.numeric(stringr::str_match(x, '[\\d\\.]+'))
  
  if(timePointType == 'D'){
    timePointMonths <- n / 30.4167
    timePointDays   <- n
  } else if(timePointType == 'M'){
    timePointMonths <- n
    timePointDays   <- n * 30.4167
  } else if(timePointType == 'Y'){
    timePointMonths <- n * 12
    timePointDays   <- n * 365
  } else {
    message('Warning - could not determine date unit for: ', x$timePointType[1])
    timePointMonths <- n
    timePointDays   <- n 
  }
  
  if(interval == 'days'){
    return(timePointDays)
  } else if(interval == 'months'){
    return(timePointMonths)
  } else {
    stop('timePoint2numeric interval error')
  }
}


termEnrichmentPlot <- function(d, n, t, file, dir = 'both'){
  d <- dplyr::arrange(d, pvalue_fdr)
  if(dir == 'down') d <- dplyr::arrange(d, genesUP - genesDOWN, pvalue_fdr)
  if(dir == 'up')   d <- dplyr::arrange(d, genesDOWN - genesUP, pvalue_fdr)
  
  d <- d[1:n,]
  
  message('min padj: ', min(d$pvalue_fdr), ' max padj: ', max(d$pvalue_fdr))
  
  # d$label <- paste0(d$term_description, '\n', signif(d$pvalue_fdr, digits=2))
  d$label <- d$term_description
  d <- dplyr::select(d, label, genesUP, genesDOWN) 
  d$genesDOWN <- d$genesDOWN * -1
  labels <- rev(unique(d$label))
  d <- tidyr::gather(d, key='var', value='val', genesUP, genesDOWN)
  d$var <- ifelse(d$var == 'genesUP', 'Increased', 'Decreased')
  d$label <- factor(d$label, levels = labels)

  p <- 
    ggplot(d, aes(label, val, fill = var)) + 
    theme_bw() +
    geom_bar(color = 'black', stat='identity') + 
    scale_fill_manual(name = 'Transcription', values = c('dodgerblue2', 'red')) +
    coord_flip() +
    labs(x = '', y = 'Genes') +
    theme(text = element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    #guides(fill=FALSE) +
    ggtitle(t)
  
  ggsave(p, file = file)
  p
}




# Function which creates a table of enriched KEGG terms from associated STRINGdb protein ids.
# The names of genes considered in the enrichments are provided and are listed with uppercase 
# letters if the gene showed increased transcription levels and listed with lowercase letters 
# if the gene showed decreased transcription levels.

createKEGGenrichmentTable <- function(o){
  e <- string_db$get_enrichment(o$STRING_id, category = 'KEGG', methodMT = "fdr", iea = FALSE)
  e <- subset(e, pvalue_fdr <= 0.05)
  
  # Retrieve all gene names for returned pathway ids.
  p <- string_db$get_term_proteins(e$term_id)
  
  # Limit pathway genes to genes in the provided data frame.
  p <- p[p$STRING_id %in% o$STRING_id,]
  
  # Create vectors of up and down regulated genes (STRINGdb ids).
  o.up   <- subset(o, log2FoldChange >= 0)$STRING_id
  o.down <- subset(o, log2FoldChange < 0)$STRING_id
  
  # Change capitalization of gene name based on up or down regulation.
  p <- dplyr::rowwise(p) %>% dplyr::mutate(preferred_name2 = ifelse(STRING_id %in% o.up, toupper(preferred_name), tolower(preferred_name))) %>% dplyr::ungroup()
  
  p <- dplyr::group_by(p, term_id) %>% 
       dplyr::summarise(genesUP   = sum(STRING_id %in% o.up),
                        genesDOWN = sum(STRING_id %in% o.down),
                        genes = paste0(sort(preferred_name2), collapse = ', ')) %>% 
       dplyr::ungroup()
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




importSalmon <- function(salmonFiles){
  library(tximport)
  
  # Concatenate salmon outputs and create DESeq2 inputs.
  d <- bind_rows(lapply(salmonFiles, function(f){
    read.table(f, sep = '\t', header = TRUE, quote = '', comment.char = '')
  }))
  
  # Create tx2gene lookup table.
  d$Name <- as.character(d$Name)
  d <- d[!duplicated(d$Name),]
  tx2gene <- data.frame(transcript_id = d$Name, 
                        gene_id       = unlist(lapply(strsplit(d$Name, '\\|'), '[[', 2)), 
                        gene_name     = unlist(lapply(strsplit(d$Name, '\\|'), '[[', 6)))
  
  tximport(files = salmonFiles, type = "salmon", tx2gene = tx2gene)
}



processSalmon <- function(txi){
  # Create metadata table.
  samples <- data.frame(genotype = factor(unlist(lapply(strsplit(colnames(txi$counts), '_'), '[[', 1))), 
                        timePt   = factor(unlist(lapply(strsplit(colnames(txi$counts), '_'), '[[', 2))),
                        Subject  = paste0('S', unlist(lapply(strsplit(colnames(txi$counts), '_'), '[[', 3))))
  samples$repGrp <- factor(paste0(as.character(samples$genotype), '_', as.character(samples$timePt)))
  rownames(samples) <- colnames(txi$counts)
  
  
  # Create DESeq2 data set.
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi, colData = samples, design = ~ repGrp)
  ddsTxi <- DESeq2::DESeq(ddsTxi)
  
  ddsTxi
}



createUCSCintSiteAbundTrack <- function(posid, abund, subject, title='intSites', outputFile='track.ucsc', visbility = 1, position=NA){
  o <- stringr::str_match_all(posid, '([^\\+^\\-]+)([+-])(\\d+)')
  d <- data.frame(chr     = unlist(lapply(o, '[[', 2)),
                  strand  = unlist(lapply(o, '[[', 3)),
                  site    = as.integer(unlist(lapply(o, '[[', 4))),
                  subject = subject,
                  abund   = as.integer(abund))
  
  d$abundBins <- cut(d$abund, 10, labels = FALSE)
  d$site2 = d$site
  d$j = '0'
  
  file.create(file = outputFile)
  if(!is.na(position)) write(paste0('browser position ', position), file=outputFile, append = TRUE)
  
  write(paste0('track name="', title, '" description="', title, '" colorByStrand="0,0,255 255,0,0" visibility=', visbility), 
        file = outputFile, 
        append = TRUE)
  
  write.table(d[,c('chr', 'site', 'site2', 'subject', 'j', 'strand')], 
              file = outputFile, 
              sep = '\t', 
              append = TRUE, 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  
  write(paste0('track type="bedGraph" name="', title, 
               ' abund bins" color=0,128,0 visiblity=full autoScale=off viewLimits=0:10 maxHeightPixels=30:30:30 visibility=', visbility), 
        file=outputFile, 
        append = TRUE)
  
  write.table(d[,c('chr', 'site', 'site2', 'abundBins')], 
              file=outputFile, 
              sep = '\t', 
              append = TRUE, 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
}



preVsPostFreqPlot <- function(sites, daysCutOff, oncoGenes, nGenesToLabel = 5, distCutOff = 10000, mustLabel = c('xxxxxx')){
  d <- data.frame(sites) %>%
       dplyr::filter(abs(nearestFeatureDist) <= distCutOff) %>%
       dplyr::mutate(preTransplantSites = n_distinct(posid[timePointDays <= daysCutOff])) %>%
       dplyr::mutate(postTransplantSites = n_distinct(posid[timePointDays > daysCutOff])) %>%
       dplyr::group_by(nearestFeature) %>%
       dplyr::summarise(preTransplant  = round( (n_distinct(posid[timePointDays <= daysCutOff]) / preTransplantSites[1]), 50),
                        postTransplant = round( (n_distinct(posid[timePointDays > daysCutOff]) / postTransplantSites[1]), 50)) %>%
       dplyr::ungroup() %>%
       dplyr::mutate(g = paste(preTransplant, postTransplant)) %>%
       dplyr::group_by(g) %>%
       dplyr::mutate(genesPerPos = n_distinct(nearestFeature)) %>%
       dplyr::ungroup() %>%
       dplyr::mutate(prePostDiff = postTransplant - preTransplant) %>%
       dplyr::mutate(oncoGene = factor(ifelse(toupper(nearestFeature) %in% toupper(oncoGenes), 'Yes', 'No'))) %>%
       dplyr::arrange(prePostDiff)

  genesToLabel <- unlist(d[c(1:nGenesToLabel,(nrow(d)-(nGenesToLabel-1)):nrow(d)),c('nearestFeature')])
  genesToLabel <- c(mustLabel, genesToLabel)
  d$geneLabel  <- ifelse(d$nearestFeature %in% genesToLabel, d$nearestFeature, '')
  
  ggplot(d, aes(preTransplant, postTransplant, color = log10(genesPerPos), shape = oncoGene)) +
    theme_bw() +
    scale_shape_manual(name = 'Oncogene', values = c(16, 15)) +
    scale_color_gradientn(name = 'log10(Data density)', colors = c("green3", "gold2", "red")) +
    geom_point(alpha = 0.4, size = 4, stroke = 0) +
    #geom_jitter(width = 0.0001, height = 0.0001) +
    geom_abline(slope=1, intercept=0, color='blue', size=0.5) +
    guides(shape = guide_legend(override.aes = list(stroke = 2))) +
    geom_text_repel(aes(label=geneLabel), color='black', size=2.5,  direction='y',box.padding=1.0, point.padding=1.25) +
    theme(legend.position="bottom") +
    guides(shape = guide_legend(title.position = "top")) +
    guides(color=FALSE) +
    labs(x = paste0('Earlier time points (<= ', daysCutOff, ' days, ', ppNum(n_distinct(subset(sites, timePointDays <= daysCutOff)$posid)), ' sites)'),
         y = paste0('Later time points (> ', daysCutOff, ' days, ', ppNum(n_distinct(subset(sites, timePointDays > daysCutOff)$posid)), ' sites)'))
}

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


termEnrichmentPlot <- function(d, n, t){
  d <- dplyr::arrange(d, pvalue_fdr)[1:n,]
  message(min(d$pvalue_fdr), ' - ', max(d$pvalue_fdr))
  # d$label <- paste0(d$term_description, '\n', signif(d$pvalue_fdr, digits=2))
  d$label <- d$term_description
  d <- dplyr::select(d, label, genesUP, genesDOWN) 
  d$genesDOWN <- d$genesDOWN * -1
  labels <- rev(unique(d$label))
  d <- tidyr::gather(d, key='var', value='val', genesUP, genesDOWN)
  d$label <- factor(d$label, levels = labels)
  ggplot(d, aes(label, val, fill = var)) + 
    theme_bw() +
    geom_bar(stat='identity') + 
    scale_fill_manual(values = c('dodgerblue2', 'red')) +
    coord_flip() +
    labs(x = '', y = 'Genes') +
    theme(text = element_text(size=14), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill=FALSE) +
    ggtitle(t)
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



nearestGenomicFeature <- function(query, genome='hg38', side='either', geneList=NULL){
  
  if(tolower(genome) == 'hg38'){
    subject       <- gt23::hg38.refSeqGenesGRanges
    subject.exons <- gt23::hg38.refSeqGenesGRanges.exons
  } else if (tolower(genome) == 'mm9') {
    subject       <- gt23::mm9.refSeqGenesGRanges
    subject.exons <- gt23::mm9.refSeqGenesGRanges.exons
  } else if (tolower(genome) == 'susscr3') {
    subject       <- gt23::susScr3.refSeqGenesGRanges
    subject.exons <- gt23::susScr3.refSeqGenesGRanges.exons
  } else if (tolower(genome) == 'macfas5') {
    subject       <- gt23::macFas5.refSeqGenesGRanges
    subject.exons <- gt23::macFas5.refSeqGenesGRanges.exons  
  } else {
    stop('There is not refSeq table for the requested genome.')
  }
  
  if(! is.null(geneList)){
    subject <- GenomicRanges::subset(subject, toupper(name2) %in% toupper(geneList))
  }
  
  # If side is not set to either, collapse the subject ranges to single positions
  if (side %in% c("5p", "3p", "midpoint")) {
    options(warn=-1)
    if (side == "5p") subject <- GenomicRanges::flank(subject, width = -1)
    if (side == "3p") subject <- GenomicRanges::flank(subject, width = -1, start = FALSE)
    if (side == "midpoint") ranges(subject) <- IRanges(mid(ranges(subject)), width = 1)
    ###subject <- subject[-GenomicRanges:::get_out_of_bound_index(subject)]
    options(warn=0)
  }
  
  options(stringsAsFactors = FALSE)
  
  query.df  <- GenomicRanges::as.data.frame(query)
  subject.df <- GenomicRanges::as.data.frame(subject)
  
  query.df$strand <- as.character(query.df$strand)
  subject.df$strand <- as.character(subject.df$strand)
  
  
  subject.exons.df <- GenomicRanges::as.data.frame(subject.exons)
  query.df$inFeature            <- FALSE
  query.df$nearestFeature       <- 'None.found'
  query.df$nearestFeatureStrand <- 'None.found'
  query.df$inFeatureExon        <- FALSE
  query.df$inFeatureSameOrt     <- FALSE
  query.df$nearestFeatureStart  <- Inf
  query.df$nearestFeatureEnd    <- Inf
  query.df$nearestFeatureDist   <- Inf
  
  o  <- suppressWarnings(GenomicRanges::nearest(query, subject, select='all', ignore.strand=TRUE))
  
  if(length(o) > 0){
    
    createCol <- function(a, b, n){
      paste0(unique(cbind(a, b))[,n], collapse=',')
    }
    
    ### browser()
    
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(
        gene   = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 1),
        strand = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 2),
        hitStart = min(subject.df[subjectHits,][['start']]),
        hitEnd   = max(subject.df[subjectHits,][['end']])) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene, strand, hitStart, hitEnd) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    query.df[a$queryHits,]$nearestFeature       <- a$gene
    query.df[a$queryHits,]$nearestFeatureStrand <- a$strand
    query.df[a$queryHits,]$nearestFeatureStart  <- a$hitStart
    query.df[a$queryHits,]$nearestFeatureEnd    <- a$hitEnd
  }
  
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject, select='all', ignore.strand=TRUE, type='any'))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    #query.df[a$queryHits,]$inFeature <- a$gene
    query.df[a$queryHits,]$inFeature <- TRUE
  }
  
  o <- suppressWarnings(GenomicRanges::distanceToNearest(query,  subject, select='all', ignore.strand=TRUE))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::top_n(-1, distance) %>%
      dplyr::ungroup() %>%
      dplyr::select(queryHits, distance) %>% 
      dplyr::distinct() %>% 
      data.frame()
    query.df[a$queryHits,]$nearestFeatureDist <- a$distance
  }
  
  query.df$nearestFeatureBoundary <- ifelse(abs(query.df$start - query.df$nearestFeatureStart) > 
                                              abs(query.df$start - query.df$nearestFeatureEnd),   
                                            query.df$nearestFeatureEnd,  
                                            query.df$nearestFeatureStart)
  
  query.df$nearestFeatureDist <- query.df$nearestFeatureDist * sign(query.df$start - query.df$nearestFeatureBoundary)
  query.df$nearestFeatureDist <- ifelse(query.df$nearestFeatureStrand=='+', query.df$nearestFeatureDist, query.df$nearestFeatureDist * -1)
  
  query.df$nearestFeatureStart    <- NULL
  query.df$nearestFeatureEnd      <- NULL
  query.df$nearestFeatureBoundary <- NULL
  
  # In exon
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject.exons, select='all', ignore.strand=TRUE, type='any'))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.exons.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    #query.df[a$queryHits,]$inFeatureExon  <- a$gene
    query.df[a$queryHits,]$inFeatureExon  <- TRUE
  }
  
  # In TU ort
  # There may be cases where a site overlaps two or more features which have the sampe orientation.
  # ie. +,+,+ and we want to reduce these down to a single unique sign for comparison.
  a <- query.df[is.na(query.df$inFeature),]  
  b <- query.df[! is.na(query.df$inFeature),]
  
  if(nrow(a) > 0 && nrow(b) > 0){
    b$nearestFeatureStrandCmp <- unlist(lapply(strsplit(b$nearestFeatureStrand, ','), function(x){ paste(unique(x), collapse=',')}))
    b$inFeatureSameOrt <- b$strand == b$nearestFeatureStrandCmp
    b$nearestFeatureStrandCmp <- NULL
    query.df <- dplyr::bind_rows(a, b)
  }
  
  GenomicRanges::makeGRangesFromDataFrame(query.df, keep.extra.columns = TRUE)
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
    geom_abline(slope=1, intercept=0, color='blue', size=0.5) +
    guides(shape = guide_legend(override.aes = list(stroke = 2))) +
    geom_text_repel(aes(label=geneLabel), color='black', size=2.5,  direction='y',box.padding=1.0, point.padding=1.25) +
    theme(legend.position="bottom") +
    guides(shape = guide_legend(title.position = "top")) +
    guides(color=FALSE) +
    labs(x = paste0('Earlier time points (<= ', daysCutOff, ' days, ', ppNum(n_distinct(subset(sites, timePointDays <= daysCutOff)$posid)), ' sites)'),
         y = paste0('Later time points (> ', daysCutOff, ' days, ', ppNum(n_distinct(subset(sites, timePointDays > daysCutOff)$posid)), ' sites)'))
}

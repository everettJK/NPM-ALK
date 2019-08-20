library(tidyverse)
library(RMySQL)
library(gt23)    # devtools::github install('everettJK/gt23')
source('./lib.R')
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE) 

salmonCommand <- '/home/everett/software/salmon-0.11.3/bin/salmon quant -p8 -i /home/everett/data/sequenceDatabases/Salmon/GRCh38.gencode29/ -l A'
runSalmonRun1 <- FALSE
runSalmonRun2 <- FALSE


# Run SALMON data set1 if requested.
#--------------------------------------------------------------------------------------------------
if(runSalmonRun1){
  if(! dir.exists('salmonOutput')) dir.create('salmonOutput')
  salmonRuns1 <- read.table('data/RNAseq_data1.tsv', header = TRUE, sep = '\t') %>%
    dplyr::rowwise() %>% 
    dplyr::select(filename, barcode) %>% 
    dplyr::mutate(command = paste0(salmonCommand, 
                                   ' -1 ', list.files('data/RNAseq_data1', pattern = barcode, full.names = TRUE)[1], 
                                   ' -2 ', list.files('data/RNAseq_data1', pattern = barcode, full.names = TRUE)[2],
                                   ' -o salmonOutput/', filename)) %>%
    dplyr::ungroup()
  invisible(sapply(salmonRuns1$command, system))
}



# Run SALMON data set2 if requested.
#--------------------------------------------------------------------------------------------------
if(runSalmonRun2){
  if(! dir.exists('salmonOutput')) dir.create('salmonOutput')
  salmonRuns2 <- read.table('data/RNAseq_data2.tsv', header = TRUE, sep = '\t') %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(command = paste0(salmonCommand, 
                                   ' -r ', dataPath,
                                   ' -o salmonOutput/', sample)) %>%
    dplyr::ungroup()
  invisible(sapply(salmonRuns2$command, system))
}



# Import Salmon RNAseq results.
#--------------------------------------------------------------------------------------------------

# Import and name salmon output files.
# The two RNAseq runs, which use different strategies, can be separated by participant numbers.

files <- list.files(path = 'salmonOutput', pattern = 'quant.sf', recursive = TRUE, full.names = TRUE)
names(files) <- unlist(lapply(strsplit(files, '\\/'), '[[', 2))

salmon1.files <- files[as.integer(str_extract(names(files), '\\d+$')) <= 3]

# Exclude all day 9 results from subject 4 per Jan's instructions.
salmon2.files <- files[as.integer(str_extract(names(files), '\\d+$')) > 3]
salmon2.files <- salmon2.files[! grepl('D9_don4', names(salmon2.files))]



# Create RNAseq data import objects.
salmon1.txi <- importSalmon(salmon1.files)
salmon2.txi <- importSalmon(salmon2.files)



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




# Import intSite data and create intSite data file.
#--------------------------------------------------------------------------------------------------

# Read in sample data.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="NPM_ALK"')
dbDisconnect(dbConn)


# Retrieve and standardize fragments, call intSites, calculate abundances, and annotate sites.
intSites <- getDBgenomicFragments(samples = samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
            stdIntSiteFragments() %>%
            collapseReplicatesCalcAbunds() %>%
            annotateIntSites()


save(list = c('tx2gene', 'salmon1.txi', 'salmon2.txi', 'intSites'), file = 'data/data.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)


# Create a gene alias data object using HGNC data in order to annotate additional STRINGdb records.

if(! file.exists('data/string_db.alt.aliases.rds')){
  library(STRINGdb)
  library(parallel)
  CPUs <- 25
  STRINGdb_dataFiles <- 'data/STRINGdb'
  
  cluster <- makeCluster(CPUs)
  HGNC <- read.table('data/HGNC.txt', sep = '\t', comment.char = '', quote = '', header = TRUE)
  string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory = STRINGdb_dataFiles)
  
  HGNC.ids <- 
    dplyr::rowwise(HGNC) %>%
    dplyr::summarise(ids = list(gsub('\\s', '', 
                                     toupper(c(Approved.Symbol, 
                                               unlist(strsplit(Previous.Symbols, ',')), 
                                               unlist(strsplit(Synonyms, ','))))))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n = ntile(1:n(), CPUs))
  
  clusterExport(cluster, varlist = c('STRINGdb_dataFiles', 'HGNC.ids'))
  
  string_db.alt.aliases <- bind_rows(parLapply(cluster, split(HGNC.ids, HGNC.ids$n), function(x){
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
  
  stopCluster(cluster)
  
  saveRDS(string_db.alt.aliases, file = 'data/string_db.alt.aliases.rds')
}

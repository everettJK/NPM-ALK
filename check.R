library(gt23)
library(dplyr)
library(GenomicRanges)

intSites <- getDBgenomicFragments('GTSP3005', 'specimen_management', 'intsites_miseq') %>%
            stdIntSiteFragments() %>%
            collapseReplicatesCalcAbunds() %>%
            annotateIntSites()
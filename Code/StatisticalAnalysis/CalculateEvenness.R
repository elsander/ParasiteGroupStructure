source('MutualInformation.R')

library(dplyr)
library(stringr)
library(vegan)
library(hash)

BuildEvennessDF <- function(pathToDataDirs){
    ## pathToDataDirs is a path to the directories containing
    ## group model partitions.
    evennessDF <- data.frame(web = character(0),
                             parcon = character(0),
                             DC = character(0),
                             G = numeric(0),
                             evenness = numeric(0))
    partHash <- hash::hash()
    ## list all directories.
    dirs <- list.dirs(pathToDataDirs, recursive = FALSE)
    for(currdir in dirs){
        ## find all partition files
        fs <- list.files(currdir, pattern = ".*Marginal.*", full.names = TRUE)
        for(f in fs){
            partition <- as.vector(as.matrix(read.table(f, header = FALSE)))
            ## Pielou's evenness
            partEvenness <- vegan::diversity(partition, index = 'shannon')/
                log(length(partition))
            partName <- stringr::str_split(basename(f), '-')[[1]]
            ## make degree-corrected flag more human-readable
            if(partName[10] == '0'){
                partDC <- 'Uncorrected'
            } else {
                partDC <- 'Corrected'
            }
            ## store evenness and partition information in a data frame
            dfLine <- data.frame(web = partName[2],
                                 parcon = partName[3],
                                 DC = partDC,
                                 G = as.numeric(partName[8]),
                                 evenness = partEvenness)
            evennessDF <- evennessDF %>% dplyr::bind_rows(dfLine)
            ## store the partition vector itself in a dictionary
            partHash[stringr::str_c(as.vector(as.matrix(dfLine[1,1:4])),
                                    collapse = '-')] <- partition
        }
    }

    ## spread data frame so that each row contains a "Corrected" and an
    ## "Uncorrected" column for evenness. This sets the data up nicely
    ## for a paired t-test.
    evennessDF2 <- tidyr::spread(evennessDF, key = DC, value = evenness)
    return(list(even = evennessDF2, dict = partHash))
}

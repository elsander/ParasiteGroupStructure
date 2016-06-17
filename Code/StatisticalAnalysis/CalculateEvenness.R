library(dplyr)
library(stringr)
library(vegan)

BuildEvennessDF <- function(pathToDataDirs){
    ## pathToDataDirs is a path to the directories containing
    ## group model partitions.
    evennessDF <- data.frame(web = character(0),
                             parcon = character(0),
                             DC = character(0),
                             G = numeric(0),
                             evenness = numeric(0))
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
            if(partName[8] == '0'){
                partDC <- 'Uncorrected'
            } else {
                partDC <- 'Corrected'
            }
            ## store evenness and partition information in a data frame
            dfLine <- data.frame(web = partName[2],
                                 parcon = partName[3],
                                 DC = partDC,
                                 G = as.numeric(partName[6]),
                                 evenness = partEvenness)
            evennessDF <- evennessDF %>% dplyr::bind_rows(dfLine)
        }
    }

    ## spread data frame so that each row contains a "Corrected" and an
    ## "Uncorrected" column for evenness. This sets the data up nicely
    ## for a paired t-test.
    evennessDF2 <- tidyr::spread(evennessDF, key = DC, value = evenness)
    return(evennessDF2)
}

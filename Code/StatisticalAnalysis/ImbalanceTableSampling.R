library(tidyr)
library(dplyr)
library(stringr)

ImbalanceTableSampling <- function(resultdir){
    fs <- list.files(resultdir)
    nfs <- length(fs)
    imbalanceFull <- data.frame(classes = character(0),
                                imbalance = numeric(0),
                                probimbal = numeric(0),
                                sigstars = character(0),
                                web = character(0),
                                parcon = character(0),
                                DC = integer(0),
                                numgroups = integer(0))

    for(i in 1:nfs){
        f <- fs[i]
        splitName <- stringr::str_split(f, '-')[[1]]
        ## if the file name does not have the correct number
        ## of pieces of information, skip it
        if(length(splitName) != 5) next
        web <- splitName[1]
        parcon <- splitName[2]
        DC <- splitName[3]
        ngrps <- splitName[4]
        imbclass <- splitName[5]
        idata <- system(stringr::str_c('tail -n 1 ', file.path(resultdir, f)),
                        intern = TRUE) %>%
                            stringr::str_trim(side = 'both') %>%
                                stringr::str_split(',')
        idata <- idata[[1]]
        ## if the file has the wrong number of columns, skip it
        if(length(idata) != 4) next
        imb <- idata[2]
        probImb <- idata[3]
        tmp <- data.frame(classes = imbclass,
                          imbalance = as.numeric(imb),
                          probimbal = as.numeric(probImb),
                          web = web,
                          parcon = parcon,
                          DC = as.integer(DC),
                          numgroups = as.integer(ngrps),
                          stringsAsFactors = FALSE)
        imbalanceFull <- dplyr::bind_rows(imbalanceFull, tmp)
    }
    imbalanceFull <- dplyr::rename(imbalanceFull, Imbalance = imbalance)
    imbalanceFull$parcon[imbalanceFull$parcon == 'Par'] <- 'No_Concomitant'
    imbalanceFull$parcon[imbalanceFull$parcon == 'ParCon'] <- 'Concomitant'
    imbalanceFull$DC <- as.character(imbalanceFull$DC)
    imbalanceFull$DC[imbalanceFull$DC == '0'] <- "Uncorrected"
    imbalanceFull$DC[imbalanceFull$DC == '1'] <- "Corrected"
    
    return(imbalanceFull)
}

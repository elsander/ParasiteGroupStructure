library(xtable)
library(tidyr)
library(dplyr)
library(stringr)

GetSigstar <- function(probimbal){
    probimbal <- as.numeric(probimbal)
    if(probimbal < .001) return("^{***}")
    if(probimbal < .01) return("^{**}")
    if(probimbal < .05) return("^{*}")
    ## if those conditions aren't met, it's not significant
    return("")
}

ImbalanceTableSampling <- function(resultdir, tabledir){
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
        ## if the file has the wrong number of columns, skip it
        if(length(idata) != 4) next
        idata <- idata[[1]]
        imb <- idata[2]
        probImb <- idata[3]
        tmp <- data.frame(classes = imbclass,
                          imbalance = as.numeric(imb),
                          probimbal = as.numeric(probImb),
                          sigstars = GetSigstar(probImb),
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

DCttest <- function(){
    out <- ImbalanceTableSampling('../../Results/Dunne-Final-Imbalance-Sampling',
                                  '.', makeTable = FALSE)
    out <- dplyr::rename(out, Imbalance = imbalance)
    out$parcon[out$parcon == 'Par'] <- 'No Concomitant'
    out$parcon[out$parcon == 'ParCon'] <- 'Concomitant'
    out$classes[out$classes == 'all'] <- 'All'
    out$classes[out$classes == 'Primary_Producers'] <- 'Producers'
    out$web <- as.factor(out$web)
    out$classes <- as.factor(out$classes)
    out$parcon <- as.factor(out$parcon)
    out$DC <- as.character(out$DC)
    out$DC[out$DC == '0'] <- "Uncorrected"
    out$DC[out$DC == '1'] <- "Corrected"
    out2 <- tidyr::spread(out, key = DC, value = Imbalance)
    ttestResults <- t.test(out2$Uncorrected, out2$Corrected,
                           alternative = "two.sided", paired = TRUE)
    return(ttestResults)
}

DCttestPar <- function(){
    out <- ImbalanceTableSampling('../../Results/Dunne-Final-Imbalance-Sampling',
                                  '.', makeTable = FALSE)
    out <- dplyr::rename(out, Imbalance = imbalance)
    out$parcon[out$parcon == 'Par'] <- 'No Concomitant'
    out$parcon[out$parcon == 'ParCon'] <- 'Concomitant'
    out$classes[out$classes == 'all'] <- 'All'
    out$classes[out$classes == 'Primary_Producers'] <- 'Producers'
    out$web <- as.factor(out$web)
    out$classes <- as.factor(out$classes)
    out$parcon <- as.factor(out$parcon)
    out$DC <- as.character(out$DC)
    out$DC[out$DC == '0'] <- "Uncorrected"
    out$DC[out$DC == '1'] <- "Corrected"
    out2 <- tidyr::spread(out, key = DC, value = Imbalance) %>%
        dplyr::filter(classes %in% c('Parasites', 'Other'))
    out2Par <- out2 %>%
        dplyr::filter(parcon == 'No Concomitant')
    out2ParCon <- out2 %>%
        dplyr::filter(parcon == 'Concomitant')
    ttestPar <- t.test(out2Par$Uncorrected, out2Par$Corrected,
                       alternative = "two.sided", paired = TRUE)
    ttestParCon <- t.test(out2ParCon$Uncorrected, out2ParCon$Corrected,
                          alternative = "two.sided", paired = TRUE)
    return(list(par = ttestPar, parcon = ttestParCon))
}

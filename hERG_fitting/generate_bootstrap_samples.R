# File:         generate_bootstrap_samples.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to generate bootstrap samples from hERG channel
#               fractional block data for 12 CiPA training drugs.
#               New drugs can be specified by running this script with the
#               command line option "-d DRUG".
#               The required input file(s) are comma-separated value (CSV)
#               files located at:
#                   data/DRUG.csv
#               For help with other options, run this script with command
#               line option "-h".
#

#--- load libraries
library(boot)
library(vroom)
library(dplyr, quietly = T, warn.conflicts = F)
library(stringr)

generate_file_names <- function(drugnames, dir) {
    purrr::map_chr(drugnames, function(d){
        file.path(dir, str_c(d, ".csv.xz"))
    }) %>% setNames(drugnames)
}

generate_bootstrap_samples <- function(drugnames, nboots = 2000, seed = 100) {

    #--- sample with replacement using boot package
    set.seed(seed, kind="L'Ecuyer-CMRG")
    
    s <- .Random.seed
    
    purrr::map(drugnames, function(drug) {
        # read in Milnes cell data
        datadf<-vroom::vroom(drug, show_col_types = F) %>%
            arrange(conc,exp,sweep,time)

        # data frame of experiments
        expdf <- datadf %>% select(conc, exp) %>% distinct()
        
        # generate bootstrap samples
        boot.out <-boot (data=expdf, statistic = function(datadf, idx){
            mean(datadf$exp[idx])
        }, R=nboots, strata=as.factor(expdf$conc))
        
        # set random seed
        s <- parallel::nextRNGStream(s)
        .Random.seed <- s
        
        boot.out
    })
}

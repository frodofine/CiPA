# File:         hERG_fitting.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to fit hERG drug binding model parameters to Milnes
#               protocol data. The drug data to be fitted must be specified
#               with "-d DRUG", where "DRUG" indicates the path to the
#               fractional block data:
#                   data/DRUG.csv
#               For help with other options, run this script with command
#               line option "-h".
#

library(deSolve)
library(cmaes)
library(vroom)
library(furrr)
library(dplyr, quietly = T, warn.conflicts = F)

source("setup_hERG_fitting.R")
source("funs/deWrapper.R")

hERG_fitting <- function(
    drug,
    drug_boots,
    boot_num,
    params = NULL,
    seed = 100,
    cores = 1,
    ctl_list=list(
        lambda = 80,
        stop.tolx = 0.001
    )
    ) {
    
    ctl_list$vectorized <- TRUE
    
    #---setup script
    Rcpp::sourceCpp("models/hergmod.cpp")
    
    print(sprintf("Running %s bootstrap sample %d on %d cores", drug, boot_num, cores))
    
    #--- set seed for reproducibility
    set.seed(seed, kind="L'Ecuyer-CMRG")
    
    setup <- setup_hERG_fitting(drug, boot_num, drug_boots)
    
    #--- initial parameter values
    if(!is.null(params)){
        # start from optimal values
        initpar<-encode_pars(params[setup$pnames], setup$pmax, setup$pp$Low, setup$pp$High)
    }else{
        # pnames,pmax from setupfile
        initpar<-runif(length(setup$pnames)) * setup$pmax
    }
    
    objfun_vec <- function(pop){
        nchr<-ncol(pop)
        chrpernode<-ceiling(nchr/cores)
        nnodes<-min(cores,nchr)
        print(Sys.time())
        
        # objfun defined in setupfile
        errors <- future_map_dfr(
            1:nnodes,
            deWrapper(chrpernode, pop, objfun, setup),
            .options = furrr_options(
                seed = NULL,
                globals = c(
                    "setup",
                    "objfun",
                    "pop",
                    "chrpernode",
                    "run_sims",
                    "decode_pars",
                    "derivs"
                )
            )
        )
        
        imin<-which.min(errors$errors)

        print(t(decode_pars(pop[,imin], setup$pmax, setup$pp$Low, setup$pp$High)))
        print(errors$errors[imin])
        
        errors$errors
    }
    
    #--- run CMA-ES
    res <- cma_es(initpar,objfun_vec,lower=0,upper=setup$pmax,control=ctl_list)
    
    #--- save best results
    pars<-signif(decode_pars(res$par, setup$pmax, setup$pp$Low, setup$pp$High), digits=4)
    names(pars) <- setup$pnames
    
    pars
}

hERG_fitting_initial <- function(
    drug,
    drug_boots,
    boot_num = 0,
    seed = 100,
    cores = 1,
    ctl_list=list(
        lambda = 80,
        stop.tolx = 0.001
    )
) {
    hERG_fitting(drug, drug_boots, boot_num, NULL, seed, cores, ctl_list)
}

hERG_fitting_boot <- function(
    drug,
    params,
    drug_boots,
    boot_num = 1:drug_boots[[drug]]$R,
    seed = 100,
    cores = 1,
    ctl_list=list(
        lambda = 80,
        stop.tolx = 0.001
    )
    ) {
    
    #--- run range of bootstraps
    purrr::map_dfr(boot_num, function(d) {
        hERG_fitting(drug, drug_boots, d, params, seed, cores, ctl_list)
    }, .id = "boot")
}

read_pars_txt <- function(file) {
    res <- read.table(file)
    
    r <- res$V2
    names(r) <- res$V1
    
    r
}

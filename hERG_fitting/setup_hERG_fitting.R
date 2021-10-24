# File:         setup_hERG_fitting.R
# Author:       Kelly Chang
#               Zhihua Li
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R helper script to setup fitting for the hERG drug binding
#               model.
#

source("funs/get_boot_data.R")
source("funs/stepprotocol.R")
source("funs/model_init.R")

encode_pars<-function(pars,pmax,low,high) pmax*log10(pars/low)/log10(high/low)
decode_pars<-function(ind,pmax,low,high) low*(high/low)^(ind/pmax)

setup_hERG_fitting <- function(drug, boot_num, drug_boots = NULL) {

    #--- code to setup simulations
    nnnfile<-"models/hergmod_drug_param_bounds.txt"
    modelname<-"hergmod"
    modeldir<-"models/"
    
    #--- load ODE model, states, and parameters
    tmp<-read.table(paste0(modeldir,modelname,"_states.txt"), col.names=c("param","value"))
    states<-setNames(tmp$value, tmp$param)
    
    tmp<-read.table(paste0(modeldir,modelname,"_pars.txt"), col.names=c("param","value"))
    pars<-setNames(tmp$value, tmp$param)
    
    #--- parameters to be fitted
    pp<-read.table(nnnfile,header=T,as.is=T)
    pnames<-pp$Parameter
    high_bounds<-pp$High
    low_bounds<-pp$Low
    
    #--- parameter encoding to 0-10 range
    pmax <- 10
    
    #--- setup model
    pars["T"] <- 37
    nbeats <- 10
    
    dummy<-stepprotocol(-80, 900, -80, 40, 0, 10000, 14060)
    fulltimes<-dummy[[1]]
    peaktimes<-dummy[[2]]
    eventdata<-dummy[[3]]
    
    hergmodel<-model_init(modelname, states, pars, pnames, fulltimes, eventdata, nbeats)
    ctlsweeps<-hergmodel$controlsweeps()
    
    #--- get bootstrap indices
    drug_boot <- drug_boots[[drug]]
    
    boot_idx<-boot::boot.array(drug_boot, indices=TRUE)
    
    if(boot_num>0){
        idx<-boot_idx[boot_num,]
    }else{
        idx<-1:ncol(boot_idx)
    }
    
    print(idx)
    
    #--- get data
    fracdata<-get_boot_data(paste0("data/",drug,".csv.xz"), idx)
    concvec<-names(fracdata) # note these are characters!
    nsweeps<-max(sapply(fracdata, length))
    if(nbeats!=nsweeps)
        stop(sprintf("Simulations have %d beats, but experiments have %d sweeps!",nbeats,nsweeps))
    
    #--- check timepoints
    haspts<-sapply(fracdata,
                   function(x){
                       sapply(1:min(length(x),nbeats),
                              function(i){
                                  all(x[[i]][,"time"]%in%ctlsweeps[[i]][,"time"])
                              })
                   })
    if(is.list(haspts)) haspts<-do.call(c,haspts) # for missing data
    if(!all(haspts)) stop("Simulation timepoints don't match data!")
    
    list(
        fracdata = fracdata,
        pnames = pnames,
        pmax = pmax,
        hergmodel = hergmodel,
        pp = pp,
        nbeats = 10
    )
}

#--- function for running drug simulations
run_sims<-function(fitpars, concvec, hergmodel){
    alldrugsweeps<-list()

    for(conc in concvec){

        drugsweeps <- hergmodel$run_drug(fitpars, as.numeric(conc)) # note conversion to numeric
        
        if(length(drugsweeps)==0)
            return(c())
        
        alldrugsweeps[[conc]]<-drugsweeps
    }
    alldrugsweeps
}

#--- define objective function
objfun<-function(ind, setup){
    if (.Platform$OS.type == "windows"){
        library(deSolve)
        source("models/hergmod.cpp.R")
    }

    # run simulations
    fracdata <- setup$fracdata
    concvec <- names(fracdata)
    
    decoded_pars <- decode_pars(ind, setup$pmax, setup$pp$Low, setup$pp$High)

    alldrugsweeps <- run_sims(decoded_pars, concvec, setup$hergmodel)

    if(length(alldrugsweeps)==0)
        return(1e50)

    ctlsweeps <- setup$hergmodel$controlsweeps()
    # compute errors
    idxV<-which(colnames(ctlsweeps[[1]])=="V")
    fval<-0
    yPred<-list()
    yO<-list()
    
    for(conc in concvec){
        yPred[[conc]]<-list()
        yO[[conc]]<-list()

        # mean squared error of picked points
        for(i in seq_along(fracdata[[conc]])){
            yO[[conc]][[i]]<-fracdata[[conc]][[i]][,"frac"]
            deptime<-fracdata[[conc]][[i]][,"time"]
            idxPeaktime<-ctlsweeps[[i]][,"time"]%in%deptime
            yPred[[conc]][[i]]<-alldrugsweeps[[conc]][[i]][idxPeaktime,"O"]/ctlsweeps[[i]][idxPeaktime,"O"]
            fval<-fval+sum((yPred[[conc]][[i]]-yO[[conc]][[i]])^2)
        }

        # negative constraints violation
        allval<-alldrugsweeps[[conc]][[setup$nbeats]][,c(-1,-idxV)]
        negativeerror<-sum(apply(allval,2,function(x) mean(pmin(x,0)^2)))
    }

    # trapping error
    numpt<- length(unlist(yO))
    significantidx<- sapply(yO, function(x) any(unlist(x)<=0.5)) # which doses have block
    firstO<-sapply(yO, function(x){sapply(x, function(z) z[1])})
    firstP<-sapply(yPred, function(x){sapply(x, function(z) z[1])})

    # handle missing data
    if(is.list(firstO)){
        tmp<-firstO
        firstO<-matrix(NA, nrow=max(sapply(tmp,length)), ncol=length(tmp))
        for(tmpi in seq_along(tmp)){
            tmp1<-tmp[[tmpi]]
            firstO[1:length(tmp1),tmpi]<-tmp1
        }
    }
    if(is.list(firstP)){
        tmp<-firstP
        firstP<-matrix(NA, nrow=max(sapply(tmp,length)), ncol=length(tmp))
        for(tmpi in seq_along(tmp)){
            tmp1<-tmp[[tmpi]]
            firstP[1:length(tmp1),tmpi]<-tmp1
        }
    }

    specerr<-0
    if(nrow(firstO)==setup$nbeats){
        yPdiff1<-(firstP[1,]-firstP[setup$nbeats,])/firstP[1,]
        yOdiff1<-(firstO[1,]-firstO[setup$nbeats,])/firstO[1,]
        yOdiff1<-pmax(0, yOdiff1)
        nabeat<-is.na(yOdiff1) | is.na(yPdiff1)
        specerr<-specerr+sum(0.2*numpt/length(concvec)*(yPdiff1-yOdiff1)[!nabeat & significantidx]^2)
    }
    if(nrow(firstO)>1){
        yPdiff2<-(firstP[1,]-firstP[2,])/firstP[1,]
        yOdiff2<-(firstO[1,]-firstO[2,])/firstO[1,]
        yOdiff2<-pmax(0,yOdiff2)
        nabeat<-is.na(yOdiff2) | is.na(yPdiff2)
        specerr<-specerr+sum(0.2*numpt/length(concvec)*(yPdiff2-yOdiff2)[!nabeat & significantidx]^2)
    }

    fval+negativeerror+specerr
}

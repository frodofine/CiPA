# File:         IC50_mcmc.R
# Author:       Kelly C. Chang
#                Zhihua Li
#               based on code by Jose Vicente
# Date:         Oct 2017
# Version:      2.0
# 
# Description:  R script to perform Markov-chain Monte Carlo (MCMC)
#               simulation for Hill equations modeling drug block of
#               ionic currents. The drug to be fitted must be specified with
#               "-d DRUG". The data for the specified drug should be located
#               in "data/drug_block.csv" (default), or at the file path
#               specified by "-f FILEPATH".
#               For help with other options, run this script with command
#               line option "-h".
#

#--- load libraries
library(FME, quietly = T)
library(coda, quietly = T)
library(dplyr, quietly = T, warn.conflicts = F)

Hill_fitting <- function(drug,
                         drug_block,
                         seednum = 100,
                         nsamp = 2000,
                         burnin = 10000,
                         thin = 10,
                         desiredChannels = "all"
  ) {
  
  sigdigits=4
  addburn<-burnin # increase burnin by addburn if convergence test fails
  
  # parameter bounds for log(IC50) and Hill coefficient (note IC50 in nM)
  lowBds<-c(log(1e-10)*0.99, 0.5*1.01)
  uppBds<-c(log(1e10)*0.99, 2.0*0.99)
  
  #full bounds include std of the error term
  fulllowBds<-c(log(1e-10), 0.5, 0)
  fulluppBds<-c(log(1e10), 2.0, 50)
  
  # Hill equation residuals
  HillRes<-function(params, datadf){
    datadf$block - Hillfun(params, datadf$conc)$block
  }
  
  # Hill equation
  Hillfun<-function(params, conc){
    data.frame(conc=conc, block=100*(1-1/(1+(conc/exp(params[1]))^params[2])))
  }
  
  # likelihood fun
  Hilllikelihood<-function(params, datadf){
    # here params has three variables: ic50, h, and std
    res<- HillRes(params[1:2], datadf)
    N<- dim(datadf)[1];  delta<- params[3] 
    -2*(-N*log(delta)-1/(2*delta^2)*sum(res^2))
    
  }
  
  SSHill<-selfStart(
    ~ 100*(1-1/(1+(conc/exp(logIC50))^h)),
    function(mCall, data, LHS) {
      xy <- sortedXyData(mCall[["conc"]], LHS, data)
      if(nrow(xy) < 4)
        stop("Too few distinct x values to fit Hill equation")
      pars<-c(logIC50=NA, h=NA)
      eps<-0.01
      if(all(xy[["y"]]<eps)){
        # 0% block for all conc
        pars[["logIC50"]]<-log(max(xy[["x"]]))
        pars[["h"]]<-1
      }else if(all(xy[["y"]]>100-eps)){
        # 100% block for all conc, should handle with care
        pars[["logIC50"]]<-log(min(xy[["x"]]))
        pars[["h"]]<-1
      }else{
        # linear fit
        ignore<- xy[["y"]]<eps | xy[["y"]]>100-eps
        xx<-log(xy[["x"]][!ignore])
        yy<-1/(1-xy[["y"]][!ignore]/100)-1
        yy<-log(yy)
        xxyy<-data.frame(xx=xx, yy=yy)
        lm.out<-lm(yy ~ xx, data=xxyy)
        cf<-coef(lm.out)
        pars[["logIC50"]]<- -cf[[1]]/cf[[2]]
        pars[["h"]]<-cf[[2]]
      }
      return(pars)
    },
    c("logIC50","h")
  )
  
  # read in patch clamp data
  datadf <- drug_block %>%
    filter(drug == !!drug) %>%
    arrange(drug, conc, channel)

  channels <- datadf$channel %>% unique %>% sort

  # save plots
  #figdir<-"figs/"
  #dir.create("figs", showWarnings = F)
  #pdf(paste0(figdir,drug,"_nls_mcmc.pdf"),width=6,height=4,useDingbats=FALSE)
  
  # save samples to table
  drugdir<-file.path("results_new",drug)
  dir.create(drugdir, recursive = T, showWarnings = F)

  if(desiredChannels == "all") {
    desiredChannels <- channels
  }

  fittedChannels <- purrr::map(desiredChannels, function(channel) {
    
    tstr<-sprintf("%s, %s channel",drug,channel)
    print(tstr)
    
    # get data frame of cells for this experiment
    expdf<-datadf %>%
      filter(drug == !!drug) %>% 
      filter(channel == !!channel) %>%
      filter(!is.na(block))
    
    # fit model using modFit (nls.lm)
    
    meanconc<-mean(expdf$conc)
    
    ic50seeds<-c(
      meanconc,meanconc*0.5,meanconc*0.3,meanconc*0.75,
      0.00001, 0.0001, 0.001, 0.01, 0.1, 1,
      100,300,500,1000,3000,5000,10000,30000,50000,100000,
      300000, 500000,1000000,3000000,5000000
    )
    
    trystarts<-cbind(log(ic50seeds),rep(0.9,length(ic50seeds)))
    
    colnames(trystarts)<-c("logIC50","h")
    tryi<-1
    
    while(tryi==1 || (inherits(tryout, "try-error") && tryi<=nrow(trystarts)) ){
      tryout<-try({
        # optim is more robust; always use it first to narrow it down unless ...
        if(!any(is.na(trystarts[tryi,]))){
          realinitpars <- optim(
            trystarts[tryi,],
            function(p) sum(HillRes(p,expdf)^2),
            lower=lowBds,
            upper=uppBds,
            method="L-BFGS-B"
          )$par
        }else{
          # if init is NA, still need to go through modFit to get tryout
          # but change init to a random pair
          realinitpars<-c(logIC50=1,h=1)
        }
        
        mf<-modFit(f=HillRes, p=realinitpars, datadf=expdf, 
                   lower=lowBds, upper=uppBds, method="Marq")
        
      })
      tryi<-tryi+1
    }
    
    # plot data only or do MCMC
    xvals<-range(expdf$conc)
    if(inherits(tryout, "try-error")){
      print("Fitting error! Skipping MCMC...")
      
      # plot data
      par(mfrow=c(1,1))
      plot(
        expdf$conc,
        expdf$block,
        log="x",  
        main=tstr,
        xlab="Concentration (nM)",
        ylab="Block (%)", 
        xlim=range(xvals),
        ylim=c(0,100)
      )
      return(NULL)
    }
    
    # variables used to plot dose-response curve
    plot_IC50<-exp(mf$par[1])
    plot_h<-mf$par[2]
    xvals<-c(xvals, plot_IC50/exp(2), plot_IC50*exp(2))
    
    # because sometimes around IC50 is not covering all concentrations tested,
    # need to expand more
    xvals<-c(xvals, min(expdf$conc)/10, max(expdf$conc)*10)
    plot_x<-exp(seq(from=min(log(xvals)), to=max(log(xvals)), length.out=100))
    plot_y<-SSHill(plot_x, log(plot_IC50), plot_h)
    
    # initialize with modFit results, as recommended in FME documentation
    startp<-setNames(mf$par, c("logIC50","h"))
    startp[3]<- sd(mf$residuals)
    
    # this will assume error variance is 1: so MCMC function won't really use
    # error term in calculating likelihood; but I'll do it myself in Hill-likelihood
  
    Covar <- NULL  
    s2prior <- NULL
    weightvar<- 0
    
    # set seed for reproducibility
    set.seed(seednum)
    burnin<-addburn

    repeat {
      
      # run MCMC
      tryout<-try({
        MCMC<-modMCMC(
          f=Hilllikelihood,
          p=startp,
          datadf=expdf, 
          lower=fulllowBds,
          upper=fulluppBds,
          jump=Covar,
          var0=s2prior,
          wvar0=weightvar,
          burninlength=burnin,
          niter=burnin+nsamp*thin,
          outputlength=nsamp,
          updatecov=10,
          ntrydr=2
        )
      })
      
      if(inherits(tryout, "try-error") || MCMC$naccepted<nsamp){
        #MCMC run still failed
        print("Error running MCMC! Removing all IC50 values...")
        return(NULL)
      }
      
      # check for convergence
      gd<-geweke.diag(as.mcmc(MCMC$pars))
      
      if(!all(abs(gd[["z"]]) < qnorm(0.975))) {
        print("Geweke diagnostic indicates lack of convergence, increasing burnin...")
        set.seed(seednum)
        burnin<-burnin+addburn
        next
      }
      
      # get sensitivity
      sR<-sensRange(
        func=Hillfun,
        parms=mf$par,
        parInput=MCMC$pars[,1:2],
        conc=plot_x,
        num=nsamp
      )

      # save results to table
      # transform back to IC50
      samples <- tibble(
        IC50 = signif(exp(MCMC$pars[,1]), digits=sigdigits),
        h = signif(MCMC$pars[,2], digits=sigdigits)
      )
      
      optimal <- tibble(
        IC50 = signif(plot_IC50, digits=sigdigits),
        h = signif(plot_h, digits=sigdigits)
      )
      
      return(
        list(
          samples = samples,
          mf = mf,
          optimal = optimal,
          MCMC = MCMC,
          sR = sR,
          plot_x = plot_x,
          plot_y = plot_y,
          expdf = expdf
        )
      )
    
    }# while not converged
  })# for channel
  
  names(fittedChannels) <- desiredChannels
  
  fittedChannels
}

plot_sensitivity <- function(channel_data, title) {
  plot(
    summary(channel_data$sR),
    quant=TRUE,
    obs=channel_data$expdf %>% select(conc,block) %>% as.matrix,
    log="x",
    xlab="Concentration (nM)",
    ylab="Block (%)",
    ylim=c(0,100),
    main=title
  )
  
  # plot modFit
  lines(channel_data$plot_x, channel_data$plot_y, col="red")
}

plot_MCMC <- function(channel_data, title) {
  plot(channel_data$MCMC$pars, main=title)
  plot(channel_data$MCMC)
  hist(channel_data$MCMC)
}

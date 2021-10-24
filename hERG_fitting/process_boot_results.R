# File:         process_boot_results.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to do postprocessing of hERG fitting results for
#               the 12 CiPA training drugs. This script aggregates
#               bootstrapping results, plots the bootstrap distributions of
#               the fitted parameters, and computes confidence intervals.
#               New drugs can be specified by running this script with the
#               command line option "-d DRUG".
#               For help with other options, run this script with command
#               line option "-h".
#

#--- load libraries
library(boot)

process_boot_results <- function(drug, drug_boots, parStart, boot_result, alpha = 0.05, sigdigits = 4) {

    #--- read in bootstraps
    boot.out <- drug_boots[[drug]]
    nboots<-boot.out$R
    
    df <- boot_result %>% select(-boot)

    #--- add slope parameter
    parStart[["slope"]]<-signif(parStart[["Kmax"]]/parStart[["halfmax"]], digits=sigdigits)
    df$slope<-signif(df$Kmax/df$halfmax, digits=sigdigits)

    #--- save samples to file
    #write.csv(df, paste0(outdir,"boot_pars.csv"), row.names=F, quote=F)

    #--- construct boot object
    boot.out$t0<-parStart
    boot.out$t<-as.matrix(df)

    #--- log transform
    islog <- !names(parStart) %in% c("n","Vhalf")
    names(parStart)[islog]<-sapply(names(parStart)[islog], function(param) paste0("log10",param))
    parStart[islog]<-log10(parStart[islog])
    df[,islog]<-log10(df[,islog])

    #--- plot and calculate quantiles
    figdir<-"figs/"
    system(paste0("mkdir -p ",figdir))
    
    pdf(paste0(figdir,drug,"_boot.pdf"),width=8,height=4.5)
    cidf<-data.frame()
    for(i in 1:length(parStart)){
        param<-names(parStart)[i]

        plot(boot.out, index=i)
        title(main=param)

        out<-quantile(df[,i], probs=c(alpha/2,1-alpha/2), na.rm=TRUE)
        names(out)<-c("lower","upper")
        out<-data.frame(param=param, value=parStart[[i]], t(out))
        cidf<-rbind(cidf,out)
    }
    dev.off()
    
    cidf
}

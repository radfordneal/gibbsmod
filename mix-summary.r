# PLOT DATA SUMMARY FROM MIXTURE MODEL RUNS.

source("plot.r")
source("mix.r")

args <- commandArgs(trailingOnly=TRUE)   # Run type and number from after --args
if (length(args)<1)
{ stop("Need run type as argument")
}
rtype <- args[1]

if (length(args)<2)
{ stop("Need number of runs as argument")
}
nruns <- as.numeric(args[2])             # Number of runs to use
runs <- 0:(nruns-1)


# PRINT SELF-TRANSITION DATA.

print (read_self (paste0("runs/mix-res-",rtype), runs))


# COMPUTE AND PLOT THE ASYMPTOTIC VARIANCES.

width <- 7
height <- 9

asv <- NULL

for (runn in runs)
{
  load (paste0("runs/mix-res-",rtype,runn))
  methods <- rownames(res)
  scans <- colnames(res)

  asymvar <- list()

  pdf(paste0("runs/mix-cl1-",rtype,runn,".pdf"),width=width,height=height)
  
  asymvar <- c (asymvar, cl1=list(
   acv_plots ("Obs 1 in cluster 1", methods, scans,
              maxlag=200.5*n, maxvar=0.12, mean_to_use=1/X,
              data=res, sub="cl1"
            )))
  
  dev.off()
  
  pdf(paste0("runs/mix-eq10-",rtype,runn,".pdf"),width=width,height=height)
  
  asymvar <- c (asymvar, eq10=list(
   acv_plots ("Obs 10 cluster size", methods, scans,
              maxlag=15.5*n, maxvar=3.5, 
              data=res, sub="eq10"
             )))
  
  dev.off()
  
  pdf(paste0("runs/mix-eq30-",rtype,runn,".pdf"),width=width,height=height)
  
  asymvar <- c (asymvar, eq30=list(
   acv_plots ("Obs 30 cluster size", methods, scans,
              maxlag=10.5*n, maxvar=7, 
              data=res, sub="eq30"
             )))
  
  dev.off()
  
  pdf(paste0("runs/mix-cl1t-",rtype,runn,".pdf"),width=width,height=height)
  
  asymvar <- c (asymvar, cl1t=list(
   acv_plots ("Obs 1 in cluster 1", methods, scans,
              maxlag=200, maxvar=0.12, mean_to_use=1/X,
              data=res, sub="cl1", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/mix-eq10t-",rtype,runn,".pdf"),width=width,height=height)
  
  asymvar <- c (asymvar, eq10t=list(
   acv_plots ("Obs 10 cluster size", methods, scans,
              maxlag=15, maxvar=3.5, 
              data=res, sub="eq10", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/mix-eq30t-",rtype,runn,".pdf"),width=width,height=height)
  
  asymvar <- c (asymvar, eq30t=list(
   acv_plots ("Obs 30 cluster size", methods, scans,
              maxlag=10, maxvar=7, 
              data=res, sub="eq30", thin=n
             )))
  
  dev.off()

  if (is.null(asv))
  { asv <- rep (list(array(0,c(dim(asymvar[[1]]),length(runs)))),
                length(asymvar))
    names(asv) <- names(asymvar)
    for (i in 1:length(asymvar))
    { rownames(asv[[i]]) <- rownames(asymvar[[i]])
      colnames(asv[[i]]) <- colnames(asymvar[[i]])
    }
  }

  for (i in 1:length(asymvar))
  { asv[[i]][,,runn+1] <- asymvar[[i]]
  }
}

print(asv)


# PLOT SUMMARIES.

pdf(paste0("runs/mix-summary-",rtype,".pdf"),width=6.5,height=8.5)

summary_plot (asv,
              methods=rownames(res),
              orders=colnames(res),
              ranges=list (cl1=c(50,300,50),        # Asymp var range & step
                           eq10=c(50,400,50),
                           eq30=c(100,500,50)),
              qnames=c(cl1="Obs 1 in cluster 1",       # Names for display
                       eq10="Obs 10 cluster size",
                       eq30="Obs 30 cluster size"))

dev.off()

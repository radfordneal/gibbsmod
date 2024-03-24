# PLOT DATA SUMMARY FROM 8x8 POTTS RUNS.

source("plot.r")
source("potts-8x8.r")

args <- commandArgs(trailingOnly=TRUE)   # Run type from after --args
if (length(args)==0)
{ stop("Need run type as argument")
}
rtype <- args[1]

if (length(args)<2)
{ stop("Need number of runs as argument")
}
nruns <- as.numeric(args[2])             # Number of runs to use
runs <- 0:(nruns-1)


# PRINT SELF-TRANSITION DATA.

print (read_self (paste0("runs/potts-res-8x8-",rtype), runs))


# COMPUTE AND PLOT THE ASYMPTOTIC VARIANCES.

width <- 7
hght <- 9

asv <- NULL

for (runn in runs)
{
  load (paste0("runs/potts-res-8x8-",rtype,runn))
  methods <- rownames(res)
  scans <- colnames(res)

  asymvar <- list()

  pdf(paste0("runs/potts-8x8-cnt1-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, cnt1=list(
   acv_plots ("Count of 1s", methods, scans,
              maxlag=32.5*n, maxvar=70, mean_to_use=n/NV,
              data=res, sub="cnt1"
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-8x8-sumsqcnt-",rtype,runn,".pdf"),
              width=width,height=hght)
  
  asymvar <- c (asymvar, sumsqcnt=list(
   acv_plots ("sum squared counts", methods, scans,
              maxlag=16.5*n, maxvar=62000, 
              data=res, sub="sumsqcnt"
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-8x8-eq-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, eq=list(
   acv_plots ("Equal neighbors", methods, scans,
              maxlag=13.5*n, maxvar=70, 
              data=res, sub="eq"
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-8x8-cnt1t-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, cnt1=list(
   acv_plots ("Count of 1s", methods, scans,
              maxlag=32, maxvar=70, mean_to_use=n/NV,
              data=res, sub="cnt1", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-8x8-sumsqcntt-",rtype,runn,".pdf"),
              width=width,height=hght)
  
  asymvar <- c (asymvar, sumsqcnt=list(
   acv_plots ("sum squared counts", methods, scans,
              maxlag=16, maxvar=62000, 
              data=res, sub="sumsqcnt", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-8x8-eqt-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, eq=list(
   acv_plots ("Equal neighbors", methods, scans,
              maxlag=13, maxvar=70, 
              data=res, sub="eq", thin=n
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


# PLOT SUMMARIES.

pdf(paste0("runs/potts-summary-8x8-",rtype,".pdf"),width=6.5,height=8.5)

summary_plot (asv,
              methods=rownames(res),
              orders=colnames(res),
              ranges=list (cnt1=c(9000,72000,10000),  # Asymp var range & step
                           sumsqcnt=c(4.5e6,42e6,10e6),
                           eq=c(5000,26000,5000)),
              qnames=c(cnt1="Count of 1s",             # Names for display
                       sumsqcnt="Sum squared counts",
                       eq="Equal neighbors"))

dev.off()

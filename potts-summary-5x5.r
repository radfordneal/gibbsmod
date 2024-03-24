# PLOT DATA SUMMARY FROM 5x5 POTTS RUNS.

source("plot.r")
source("potts-5x5.r")

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

print (read_self (paste0("runs/potts-res-5x5-",rtype), runs))


# COMPUTE AND PLOT THE ASYMPTOTIC VARIANCES.

width <- 7
hght <- 9

asv <- NULL

for (runn in runs)
{
  load (paste0("runs/potts-res-5x5-",rtype,runn))
  methods <- rownames(res)
  scans <- colnames(res)

  asymvar <- list()

  pdf(paste0("runs/potts-5x5-cnt1-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, cnt1=list(
   acv_plots ("Count of 1s", methods, scans,
              maxlag=10.5*n, maxvar=3.4, minacov=-0.5*3.4, mean_to_use=n/NV,
              data=res, sub="cnt1", rnd=2
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-5x5-sumsqcnt-",rtype,runn,".pdf"),
              width=width,height=hght)
  
  asymvar <- c (asymvar, sumsqcnt=list(
   acv_plots ("sum squared counts", methods, scans,
              maxlag=16.5*n, maxvar=120, 
              data=res, sub="sumsqcnt"
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-5x5-eq-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, eq=list(
   acv_plots ("Equal neighbors", methods, scans,
              maxlag=5.5*n, maxvar=7.8, 
              data=res, sub="eq", rnd=2
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-5x5-cnt1t-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, cnt1=list(
   acv_plots ("Count of 1s", methods, scans,
              maxlag=10, maxvar=3.4, minacov=-0.5*3.4, mean_to_use=n/NV,
              data=res, sub="cnt1", thin=n, rnd=2
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-5x5-sumsqcntt-",rtype,runn,".pdf"),
              width=width,height=hght)
  
  asymvar <- c (asymvar, sumsqcnt=list(
   acv_plots ("sum squared counts", methods, scans,
              maxlag=16, maxvar=120, 
              data=res, sub="sumsqcnt", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/potts-5x5-eqt-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, eq=list(
   acv_plots ("Equal neighbors", methods, scans,
              maxlag=5, maxvar=7.8, 
              data=res, sub="eq", thin=n, rnd=2
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

pdf(paste0("runs/potts-summary-5x5-",rtype,".pdf"),width=6.5,height=8.5)

summary_plot (asv,
              methods=rownames(res),
              orders=colnames(res),
              ranges=list (cnt1=c(10,145,20),          # Asymp var range & step
                           sumsqcnt=c(1000,9000,2000),
                           eq=c(50,450,100)),
              qnames=c(cnt1="Count of 1s",             # Names for display
                       sumsqcnt="Sum squared counts",
                       eq="Equal neighbors"))

dev.off()

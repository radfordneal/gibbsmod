# PLOT DATA SUMMARY FROM BELIEFNET RUNS.

source("plot.r")
source("beliefnet.r")

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

print (read_self (paste0("runs/beliefnet-res-",rtype), runs))


# COMPUTE AND PLOT THE ASYMPTOTIC VARIANCES.

width <- 7
hght <- 9

asv <- NULL

#marg0 <- marginal0(par,1)
marg1 <- marginal1(par,1)
marg2 <- marginal2(par,1)
joint <- joint02(par,1,1)

for (runn in runs)
{
  load (paste0("runs/beliefnet-res-",rtype,runn))
  methods <- rownames(res)
  scans <- colnames(res)

  asymvar <- list()

  pdf(paste0("runs/beliefnet-mar1-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, mar1=list(
   acv_plots ("Unit 1 of layer 1 is 1", methods, scans,
              maxlag=11.5*n, maxvar=0.19, mean_to_use=marg1[1], rnd=2,
              data=res, sub="mar1"
             )))
  
  dev.off()
  
  pdf(paste0("runs/beliefnet-mar2-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, mar2=list(
   acv_plots ("Unit 1 of layer 2 is 1", methods, scans,
              maxlag=11.5*n, maxvar=0.08, mean_to_use=marg2[1], rnd=2,
              data=res, sub="mar2"
             )))
  
  dev.off()
  
  pdf(paste0("runs/beliefnet-and02-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, and02=list(
   acv_plots ("Unit1 layer0 and unit1 layer2 is 1", methods, scans,
              maxlag=11.5*n, maxvar=0.055, mean_to_use=joint[1][1], rnd=2,
              data=res, sub="and02"
             )))
  
  dev.off()
  
  pdf(paste0("runs/beliefnet-mar1t-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, mar1=list(
   acv_plots ("Unit 1 of layer 1 is 1", methods, scans,
              maxlag=11, maxvar=0.19, mean_to_use=marg1[1], rnd=2,
              data=res, sub="mar1", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/beliefnet-mar2t-",rtype,runn,".pdf"),width=width,height=hght)
  
  asymvar <- c (asymvar, mar2=list(
   acv_plots ("Unit 1 of layer 2 is 1", methods, scans,
              maxlag=11, maxvar=0.08, mean_to_use=marg2[1], rnd=2,
              data=res, sub="mar2", thin=n
             )))
  
  dev.off()
  
  pdf(paste0("runs/beliefnet-and02t-",rtype,runn,".pdf"),
      width=width,height=hght)
  
  asymvar <- c (asymvar, and02=list(
   acv_plots ("Unit1 layer0 and unit1 layer2 is 1", methods, scans,
              maxlag=11, maxvar=0.055, mean_to_use=joint[1][1], rnd=2,
              data=res, sub="and02", thin=n
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

pdf(paste0("runs/beliefnet-summary-",rtype,".pdf"),width=6.5,height=8.5)

summary_plot (asv,
              methods=rownames(res),
              orders=colnames(res),
              ranges=list (mar1=c(2,11,2),        # Asymp var range & step
                           mar2=c(0.5,5,1),
                           and02=c(0.5,5,1)),
              qnames=c(mar1="Unit 1 of layer 1 is 1",       # Names for display
                       mar2="Unit 1 of layer 2 is 1",
                       and02="Unit 1 layer 0 and unit 1 layer 2 is 1"))

dev.off()

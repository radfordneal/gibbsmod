# MIXTURE MODEL EXPERIMENTAL RUNS - TAIL COMMANDS.
#
# Included by mix-runs-a.r, mix-runs-b.r, etc. after some variables are set.

source("methods.r")
source("scans.r")
source("mix.r")
source("plot.r")


# RUN PARAMETERS.

K <- 200000    # Number of iterations to simulate for each run (each
               # consisting of n updates)

# DO RUNS.

set.seed(1); shuffle_order <- sample(n)  # for shuffled sequential scan
cat("Shuffle order:",shuffle_order,"\n")

args <- commandArgs(trailingOnly=TRUE)   # Run number from after --args
runn <- if (length(args)==0) 0 else as.integer(args)

seed <- (runn+1)*(1000000+10000*which(rtype==c("a","b","c")))

set.seed(seed<-seed+100)

s0 <- sample(1:X,n,replace=TRUE) # Random initial mixture component assignments

res <- matrix(list(NULL),length(meth),length(scan))
rownames(res) <- names(meth)
colnames(res) <- names(scan)

for (mth in 1:length(meth))
{ for (scn in 1:length(scan))
  { cat (paste0("\n",names(meth)[mth],", ",names(scan)[scn],":\n"))
    set.seed(seed<-seed+100)
    print (system.time (
            res[[mth,scn]] <- mix_run (s0, K, meth[[mth]], scan[[scn]])))
    print (res[[mth,scn]]$final)
  }
}


# SAVE THE RESULTS FOR LATER PLOTTING AND ANALYSIS.

save (res, file=paste0("runs/mix-res-",rtype,runn), version=2)

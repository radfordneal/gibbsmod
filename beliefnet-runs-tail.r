# BELIEFNET EXPERIMENTAL RUNS - TAIL COMMANDS.
#
# Included by beliefnet-runs-a.r, beliefnet-runs-b.r, etc. after some
# variables are set.


# RUN PARAMETERS.

K <- 10000      # Number of iterations to simulate for each run (each
                # consisting of n updates)

# DO RUNS.

set.seed(1); shuffle_order <- sample(n)  # For shuffled sequential scan
cat("Shuffle order:",shuffle_order,"\n")

args <- commandArgs(trailingOnly=TRUE)   # Run number from after --args
runn <- if (length(args)==0) 0 else as.integer(args)

seed <- (runn+1)*(1000000+10000*which(rtype==c("a","b","c")))

set.seed(seed<-seed+100)

s0 <- c (sample(1:M0,N0,replace=TRUE),
         sample(1:M1,N1,replace=TRUE),
         sample(1:M2,N2,replace=TRUE))

res <- matrix(list(NULL),length(meth),length(scan))
rownames(res) <- names(meth)
colnames(res) <- names(scan)

for (mth in 1:length(meth))
{ for (scn in 1:length(scan))
  { cat (paste0("\n",names(meth)[mth],", ",names(scan)[scn],":\n"))
    set.seed(seed<-seed+100)
    print (system.time (
      res[[mth,scn]] <- beliefnet_run (par, s0, K, meth[[mth]], scan[[scn]])))
    print (res[[mth,scn]]$final)
  }
}


# SAVE DATA FOR LATER ANALYSIS.

save (res, file=paste0("runs/beliefnet-res-",rtype,runn), version=2)

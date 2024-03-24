# PRODUCE PLOTS OF POTTS MODEL GIBBS SAMPLING AND MODIFIED DISTRIBUTIONS.
#
# Assumes that the number of possible values for a site is four.

source("methods.r")

potts_cond_plot <- function (b)
{
  par(mar=c(2.5,2,2,0.23))
  par(mgp=c(1.20,0.35,0),tcl=-0.35)
  par(mfrow=c(6,4))
  par(lab=c(4,3,7))

  colours <- c("red","green","blue","orange","magenta")
  range <- c(-0.03,1.03)

  nlist <- list(c(1,2,3,4), c(1,1,2,3), c(1,1,2,2), c(1,1,1,2), c(1,1,1,1))
  mlist <- list("GS","UNAM","ZDNAM","HST","ST","ZFSS")

  pi <- function (neighbors)
  { p <- numeric(4)
    for (i in 1:4)
    { p[i] <- exp(b*sum(neighbors==i))
    }
    p/sum(p)
  }

  show_key = TRUE

  for (method in mlist)
  { trans <- get(paste0("trans_",method))
    for (curr in 1:4)
    { plot(1:4,rep(range,length=4),type="n",yaxs="i",xlim=c(0.8,4.2),
           xlab="", ylab="")
      title(paste(method,"from",curr))
      abline(h=c(0,0.5,1),col="lightgray",lty=1)
      if (show_key)
      { text(2.55,0.95,"Neighbors:",pos=4)
        for (i in 1:length(nlist))
        { lines(c(2.3,3),rep(0.93-0.08*i,2),col=colours[i])
          text(3.05,0.93-0.08*i,paste(nlist[[i]],collapse=" "),pos=4)
        }
        show_key <- FALSE
      }
      for (i in 1:length(nlist))
      { p <- trans(pi(nlist[[i]]),curr)
        lines(1:4,p,type="b",pch=19,col=colours[i])
      }
    }
  }
}

source("potts-8x8.r")
pdf("runs/potts-cond-8x8.pdf",width=10,height=14)
potts_cond_plot(b)

source("potts-5x5.r")
pdf("runs/potts-cond-5x5.pdf",width=10,height=14)
potts_cond_plot(b)

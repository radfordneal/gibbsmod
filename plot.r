# FUNCTIONS FOR READING DATA AND PLOTTING.


# FUNCTION FOR PLOTTING AUTOCOVARIANCES.  Also displays the asymptotic variance
# estimate in the plot title, and returns this as its value.

plot_auto_cov <- function (vals,maxlag,maxvar,minacov=-0.15*maxvar,rnd=0,
                           xlab="",name="",mean_to_use=mean(vals),thin=1)
{ valsm <- vals - mean_to_use
  a <- acf(valsm,lag.max=maxlag,type="covariance",demean=FALSE,
           ylab="", # ylab=paste("var",round(sum(valsm^2)/length(valsm),rnd)),
           ylim=c(minacov,maxvar), xlab=paste(xlab,"  "),main="")$acf
  asymv <- thin * (a[1]+2*sum(a[-1]))
  title (name, line=0.4)
  title (paste("asym var",round(asymv,rnd)," "), line=-1, adj=1)
  asymv
}


# PRODUCE MULTIPLE ACV PLOTS.  Takes a data matrix of results; returns
# a matrix of asymptotic variance estimates, matching the matrix of data.

acv_plots <- function (name, row_labels, col_labels, maxlag, maxvar, 
                       data, thin=1, sub, ...)
{
  nr <- length(row_labels)
  nc <- length(col_labels)

  a <- matrix(0,nr,nc)
  rownames(a) <- row_labels
  colnames(a) <- col_labels

  par(mfrow=c(nr,nc))
  par(mar=c(2.5,2,2,0.23))
  par(mgp=c(1.1,0.25,0))
  par(tcl=-0.15)
  par(cex.axis=0.8)
  par(cex.lab=0.8)
  par(cex.main=0.8)
  par(font.main=1)

  for (r in 1:nr)
  { for (c in 1:nc)
    { d <- data[[r,c]]
      if (!missing(sub)) d <- d[[sub]]
      if (thin!=1) d <- d[seq(thin,length(d),by=thin)]
      a[r,c] <- plot_auto_cov (d, maxlag, maxvar, 
                  name=name, xlab=paste0(col_labels[c],", ",row_labels[r]),
                  thin=thin, ...)
    }
  }

  a
}


# READ SELF-TRANSITION DATA.

read_self <- function (file_base, runs)
{
  self <- 0
  prhalf <- 0
  self_pr <- 0
  min_self <- 0
  for (r in 1:length(runs))
  { load(paste0(file_base,runs[r]))
    s <- as.numeric(lapply(res,function(x)x$self))
    dim(s)<-dim(res)
    dimnames(s)<-dimnames(res)
    self <- self + s
    s <- as.numeric(lapply(res,function(x)x$prhalf))
    dim(s)<-dim(res)
    dimnames(s)<-dimnames(res)
    prhalf <- prhalf + s
    s <- as.numeric(lapply(res,function(x)x$self_pr))
    dim(s)<-dim(res)
    dimnames(s)<-dimnames(res)
    self_pr <- self_pr + s
    s <- as.numeric(lapply(res,function(x)x$min_self))
    dim(s)<-dim(res)
    dimnames(s)<-dimnames(res)
    min_self <- min_self + s
  }

  list (self = self/length(runs),
        prhalf = prhalf/length(runs),
        self_pr = self_pr/length(runs),
        min_self = min_self/length(runs))
}


# PRODUCE SUMMARY PLOT.

ordcol <- c(`Random`="black",
            `Sequential`="red",
            `Shuffled Sequential`="orange",
            `Checkerboard`="green",
            `Random order`="blue",
            `Random order x4`="cyan")

summary_plot <- function (asv,methods,orders,ranges,qnames)
{
  cat("methods:",methods,"\n")
  par(mfcol=c(length(ranges),2))
  par(mar=c(2.5,2,2.5,0.2))
  par(mgp=c(1.1,0.25,0))
  par(tcl=-0.15)
  par(cex.axis=0.8)
  par(cex.main=0.9)

  qnames <- c(qnames,paste(qnames,"- thinned"))
  ranges <- rep(ranges,2)
  
  key <- function (rng)
  { for (j in 1:length(orders))
    { ypos <- rng[1]+(rng[2]-rng[1])*(1.035-j*0.042)
      lines(length(methods)*c(0.545,0.599),rep(ypos,2),col=ordcol[orders[j]])
      text(length(methods)*0.585, ypos, orders[j], pos=4, cex=0.75)
    }
  }
  
  for (i in 1:length(asv))
  { if (is.null(asv[[i]]))
    { plot (0, 0, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
      next
    }
    a <- asv[[i]]
    rng <- ranges[[i]]
    np <- length(a)
    plot (rep(1:length(methods),length=np), a,
          type = "n", ylim=rng[1:2], ylab="",
          xlim=c(0.8,length(methods)+0.2), xlab="", xaxt="n")
    abline(h=seq(0,rng[2],by=rng[3]),col="gray")
    points (rep(1:length(methods),length=np), a, pch=20, 
            col=rep(ordcol[orders],each=length(methods),length=np))
    for (j in 1:length(orders))
    { lines(1:length(methods),rowMeans(a[,j,,drop=FALSE]),col=ordcol[orders[j]])
    }
    title(qnames[i])
    key(ranges[[i]])
    for (j in 1:length(methods))
    { txt <- rownames(a)[j]
      mtext(txt, side=1, line=0.2, adj=0, at=j-nchar(txt)*0.06, cex=0.5)
    }
  }
}

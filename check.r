# FUNCTION TO CHECK WHETHER METHODS WORK AS INTENDED.


# CHECK THAT THE RESULTS FROM TWO ALGORITHMS MATCH.  The tolerance argument
# gives the slop for floating point comparisons.

check_same <- function (trans_a,trans_b,p,tolerance=1e-10)
{ for (i in 1:length(p))
  { if (any(abs(trans_a(p,i)-trans_b(p,i))>tolerance)) 
    { stop("transitions are different")
    }
  }
  invisible()
}


# MAKE TRANSITION MATRIX FROM ROWS.  Also prints various data /
# checks, unless check is FALSE.  Arguments are the function to make
# rows of the matrix, and the Gibbs sampling probabilities.  For use
# in debugging and testing.

mat <- function (trans, p, check=TRUE)
{ m <- length(p)
  mat <- NULL
  for (i in 1:m) mat <- rbind (mat, trans(p,i))
  if (check)
  { cat("norm:",round(sum(p),6),
        "  sum:",round(rowSums(mat),6),
        "  check:",round(p%*%mat-p,6),
        "  self:",round(sum(p*diag(mat)),6),"\n")
  }
  mat
}


# CHECK THAT TRANSITIONS ARE VALID.  Optionally also checks that they
# do not have unnecessary non-zero self transitions, and that they are
# reversible.  Can also check that they are the same as some other
# transition, or the reverse of some other transition. The tolerance
# argument gives the slop for floating point comparisons.

check <- function (trans, p, zero_self=FALSE, reversible=FALSE, 
                   same_as=NULL, reverse_of=NULL, 
                   nonpos_eigen=FALSE, neg_eigen=FALSE, 
                   eigen_small_as=NULL, cov_dom=NULL,
                   tolerance=1e-10)
{
  if (any(p<0) || any(p>1) || abs(sum(p)-1)>tolerance)
  { stop("invalid probabilities")
  }

  M <- mat(trans,p,check=FALSE)
  if (nrow(M)!=length(p) || ncol(M)!=length(p))
  { stop("invalid matrix dimensions")
  }

  if (any(M<0))
  { stop("negative transition probability")
  }

  if (any(M>1))
  { stop("transition probability greater than one")
  }

  if (any(rowSums(M)-1>tolerance))
  { stop("transitions from a state do not sum to one")
  }

  if (any(abs(p-p%*%M)>tolerance))
  { stop("transitions do not leave distribution invariant")
  }

  if (zero_self)
  { if (any(diag(M)!=0) && all(p<=0.5-tolerance))
    { stop("non-zero self transition when probability not greater than 1/2")
    }
  }

  if (reversible)
  { if (any(abs(p*M-t(p*M))>tolerance))
    { stop("transitions are not reversible")
    }
  }

  if (!is.null(same_as))
  { M2 <- mat(same_as,p,check=FALSE)
    if (any(abs(M-M2)>tolerance)) stop("transitions not the same")
  }

  if (!is.null(reverse_of))
  { M2 <- mat(reverse_of,p,check=FALSE)
    if (any(abs(p*M-t(p*M2))>tolerance)) stop("transitions not reversals")
  }

  if (nonpos_eigen || neg_eigen || !is.null(eigen_small_as))
  { v <- eigen(M,only.values=TRUE)$values[-1]
    if (any(abs(Im(v))>1e-7))
    { stop("complex eigenvalue")
    }
    v <- Re(v)
    if (nonpos_eigen && any(v>1e-7))
    { stop(paste("eigenvalue positive:",paste(round(v,7),collapse=" ")))
    }
    if (neg_eigen && any(v>-1e-7))
    { stop(paste("eigenvalue not negative:",paste(round(v,7),collapse=" ")))
    }
    if (!is.null(eigen_small_as))
    { M2 <- mat(eigen_small_as,p,check=FALSE)
      u <- eigen(M2,only.values=TRUE)$values[-1]
      if (any(abs(Im(u))>1e-7))
      { stop("complex eigenvalue for other")
      }
      u <- Re(u)
      if (any(v>u+1e-7))
      { stop("eigenvalues not smaller than other")
      }
    }
  }
  if (!is.null(cov_dom))
  { M2 <- mat(cov_dom,p,check=FALSE)
    v <- eigen (diag(p)%*%(M2-M),only.values=TRUE)$values
    if (any(abs(Im(v))>1e-7))
    { stop("complex eigenvalue")
    }
    v <- Re(v)
    if (any(v < -1e-7))
    { stop(paste("Doesn't covariance dominate other:",
                 paste(round(v,7),collapse=" ")))
    }
  }

  invisible()
}


# DO CHECKS ON RANDOM DISTRIBUTIONS.  Also does a few manual checks first.

random_checks <- 
  function(trans, zero_self=FALSE, reversible=FALSE, 
           same_as=NULL, reverse_of=NULL,
           nonpos_eigen=FALSE, neg_eigen=FALSE, 
           eigen_small_as=NULL, cov_dom=NULL,
           tolerance=1e-10, seed=1)
{
  check(trans,c(1),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)

  check(trans,c(0,1),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(1,0),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.5,0.5),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.3,0.7),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.7,0.3),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)

  check(trans,c(0,1,0),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.5,0,0.5),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.6,0.2,0.2),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.6,0.1,0.3),zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)

  check(trans,c(0,0.5,0,0.5),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0,0.5+1e-11,0,0.5+2e-11),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0,0.5-1e-11,0,0.5-1e-11),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0,0.5-1.2e-11,0,0.5-1.0e-11),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0,0.3,0,0.3,0,0.4),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.1,0.3,0.3,0.3),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.20,0.10,0.30,0.40),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.20,0.20,0.20,0.40),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.20,0.09,0.31,0.40),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.05,0.10,0.40,0.45),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.05,0.30,0.31,0.34),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.05,0.20,0.35,0.40),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.05,0.10,0.25,0.30,0.30),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.03,0.06,0.12,0.19,0.29,0.31),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.04,0.07,0.10,0.19,0.29,0.31),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.05,0.05,0.1,0.1,0.1,0.3,0.3),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
  check(trans,c(0.03,0.04,0.11,0.11,0.11,0.3,0.3),
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)

  check(trans,c(1-1e-14,2e-13,0.5e-13,1e-13),  # probs have some round-off error
        zero_self,reversible,same_as,reverse_of,
        nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)

  set.seed(seed)
  for (m in 1:7)
  { for (i in 1:1000)
    { q <- runif(m,0.0001,1)
      p <- q/sum(q)
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      p <- exp(3*q)/sum(exp(3*q))
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      p <- -log(q)/sum(-log(q))
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      if (m<3) next
      q[m] <- q[m-1]
      p <- q/sum(q)
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      p <- exp(3*q)/sum(exp(3*q))
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      p <- -log(q)/sum(-log(q))
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      if (m<4) next
      q[1] <- q[2]
      p <- q/sum(q)
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      p <- exp(3*q)/sum(exp(3*q))
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
      p <- -log(q)/sum(-log(q))
      check(trans,p,zero_self,reversible,same_as,reverse_of,
            nonpos_eigen,neg_eigen,eigen_small_as,cov_dom,tolerance)
    }
  }
}

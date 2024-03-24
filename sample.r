# METHODS FOR SAMPLING USING MODIFIED GIBBS SAMPLING TRANSITION PROBABILITIES.


# SAMPLE FROM THE DISTRIBUTION FOR THE NEXT STATE WITH MHGS.
# Uses Algorithm 3.

sample_trans_MHGS <- function (i,samp,prob,epsilon=0)
{ 
  j <- samp()

  pi <- prob(i)
  if (pi>1-epsilon)
  { return(j)
  }

  pj <- prob(j)
  if (pj>1-epsilon)
  { return(j)
  }

  if (j!=i && pj>=pi)
  { return(j)
  }

  if (i==j)
  { seen <- i
    sum <- pi
  }
  else
  { seen <- c(i,j)
    sum <- pi + pj
  }

  h <- j
  ph <- pj

  while (sum<epsilon)
  { k <- samp()
    if (! k %in% seen)
    { pk <- prob(k)
      if (pk>1-epsilon)
      { return(j)
      }
      sum <- sum + pk
      seen <- c(seen,k)
      if (h==i)
      { h <- k
        ph <- pk
      }
    }
  }

  return (if (runif(1)<(1-pi)/(1-ph)) h else i)
}


# SAMPLE FROM THE DISTRIBUTION OF THE NEXT STATE WITH ITERATED METROPOLIZED GS.

sample_trans_IMGS <- function (i,samp,prob,m)
{ 
  pi <- prob(i)
  j <- samp()
  pj <- prob(j)
  p <- NULL

  while (j==i || pj<pi)
  {
    if (is.null(p))
    { p <- numeric(m)
      for (k in 1:m) 
      { p[k] = if (k==i) pi else if (k==j) pj else prob(k)
      }
      o = order(p)
      if (i==o[m])
      { t <- trans_1(p,i,o)
        return (sample(m,1,prob=t))
      }
      t <- numeric(m+1)
      t[m+1] <- 0
      k <- m
      while (o[k]!=i)
      { t[k] <- t[k+1] + p[o[k]]
        k <- k-1
      }
      low <- i+1
      f <- numeric(m)
      f[i] <- 1
    }

    if (j<low)
    { while (k!=j)
      { t[k] <- t[k+1] + p[o[k]]
        f[k-1] = f[k] * (t[k] - p[o[k]]) / (t[k] - p[o[k-1]])
        k <- k-1
      }
      low <- j
    }
  }

  return(j)
}

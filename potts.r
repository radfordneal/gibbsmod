# POTTS MODEL.


# FIND NUMBER OF EQUAL NEIGHBORS IN A POTTS MODEL STATE.  The state, s, 
# is a 2D array of integer values (from 1).

equal_neighbors <- function (s)
{ nr <- nrow(s)
  nc <- ncol(s)
  e <- 0L
  for (i in 1:nr)
  { for (j in 1:nc)
    { v <- s[i,j]
      e <- e + (v==s[i%%nr+1,j]) + (v==s[i,j%%nc+1])
    }
  }
  e
}


# FIND THE VALUES OF NEIGHBORS.  The state, s, is a 2D array of
# integer values and (i,j) is the location of the site to look at.
# The value returned is the vector of values for the four neighbors of
# this site (above, below, left, right, wrapping around).

neighbor_values <- function (s, i, j)
{ nr <- nrow(s)
  nc <- ncol(s)
  c (s[i%%nr+1,j], s[(i+nr-2)%%nr+1,j], s[i,j%%nc+1], s[i,(j+nc-2)%%nc+1])
}


# FIND THE CONDITIONAL DISTRIBUTION FOR A SITE IN A POTTS MODEL.  The
# state, s, is a 2D array of integer values (from 1) for sites, b is
# the bonding strength, and (i,j) is the location of the site.  The
# value returned is a vector of probabilities for possible values at
# the site.  The values of neighbors are given as the final argument,
# with the default being to compute them here.

cond <- function (s, b, i, j, nv=neighbor_values(s,i,j))
{ e <- numeric(NV)
  for (v in 1:NV)
  { e[v] <- -b * sum(nv==v)
  }
  p <- exp(-e)
  p/sum(p)
}


# POTTS SAMPLING RUN.  Does K iterations, each of which consists of
# updates of sites starting with initial state 's', using function
# 'indexes' (taking number of variables as argument) to produce the
# vector of variable indexes to use in an iteration (must always be
# the same length, default sequential), and using the function
# 'method' to modify the Gibbs sampling probabilities (default, no
# modification).  The value returned is a list with vectors 'eq' being
# the numbers of equal neighbors, after each update, 'cnt1' being the
# count of sites with value 1, after each update, 'sumsqcnt' being the
# sum of the squares of counts for each value, after each update,
# 'self' being the frequency of self transitions, 'self_pr' being the
# average probabilities of a self transitions, 'min_self' being the
# average minimum probabilities of self transitions, and 'final' being
# the final state.  Prints the self transition statistics, and the
# means of 'eq' and 'cnt1'.  Note that 'eq' and counts needed for
# 'cnt1' and 'sumsqcnt' are computed at the beginning, but then
# updated incrementally after each variable update.

potts_run <- function (s, b, K, method=function(p,v)p, indexes=function(n)1:n)
{ nr <- nrow(s)
  nc <- ncol(s)
  n <- nr*nc
  ix <- indexes(n)
  eq <- integer(K*n)
  cnts <- matrix(integer(K*n*NV),NV,K*n)
  curr_eq <- equal_neighbors(s)
  curr_cnts <- integer(NV)
  for (i in 1:NV) curr_cnts[i] <- sum(s==i)
  self_trans <- 0
  self_trans_pr <- 0
  min_self_trans_pr <- 0
  prhalf <- 0
  for (t in 1:K)
  { for (k in 1:n)
    { w <- ix[k]
      i <- (w-1)%/%nc + 1
      j <- (w-1)%%nc + 1
      nv <- neighbor_values(s,i,j)
      v <- s[i,j]
      pi <- cond(s,b,i,j,nv)
      mxpi <- max(pi)
      if (mxpi>=0.5) prhalf <- prhalf+1
      msp <- if (mxpi<=0.5) 0 else 2*mxpi-1
      min_self_trans_pr <- min_self_trans_pr + msp
      p <- method(pi,v)
      self_trans_pr <- self_trans_pr + p[v]
      new <- sample(length(p),1,prob=p)
      if (new==v)
      { self_trans <- self_trans + 1
      }
      else
      { s[i,j] <- new
        curr_eq <- curr_eq - sum(nv==v) + sum(nv==new)
        curr_cnts[v] <- curr_cnts[v] - 1
        curr_cnts[new] <- curr_cnts[new] + 1
      }
      h <- (t-1)*n + k
      eq[h] <- curr_eq
      cnts[,h] <- curr_cnts
    }
    ix <- indexes(n)
  }

  cat("self transition frequency:",self_trans/(K*n),"\n")
  cat("average self transition probability:",self_trans_pr/(K*n),"\n")
  cat("minimum self transition probability:",min_self_trans_pr/(K*n),"\n")
  cat("fraction of time max pr > 1/2:",prhalf/(K*n),"\n")

  cat("\nequal neighbors: mean",mean(eq)," variance",var(eq),"\n")
  print(table(eq))
  cnt1 <- cnts[1,]
  cat("\ncount of 1 values: mean",mean(cnt1)," variance",var(cnt1),"\n\n")
  sumsqcnt <- colSums(cnts^2)
  cat("sum sq counts: mean",mean(sumsqcnt)," variance",var(sumsqcnt),"\n\n")

  cat("final counts:",curr_cnts,"\n\n")
  list (eq=eq, cnt1=cnt1, sumsqcnt=sumsqcnt, 
        self=self_trans/(K*n), prhalf=prhalf/(K*n),
        self_pr=self_trans_pr/(K*n), min_self_pr=min_self_trans_pr/(K*n),
        final=s)
}

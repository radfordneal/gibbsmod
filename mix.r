# BAYESIAN MIXTURE MODEL.


# PARAMETERS.

n <- 30    # Number of observations
M <- 10    # Number of binary variables in each observation

X <- 9     # Number of mixture components


# DATA.

data <- c(

  1, 1, 1, 1, 0, 0, 0, 0, 1, 0,
  1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0, 0, 1, 0,
  1, 0, 1, 1, 0, 0, 0, 0, 1, 0,
  1, 1, 1, 1, 0, 0, 0, 0, 0, 1,
  1, 1, 1, 1, 0, 0, 1, 0, 1, 1,
  0, 1, 1, 1, 0, 0, 0, 0, 0, 0,

  0, 0, 0, 0, 1, 1, 1, 1, 1, 0,
  0, 0, 0, 0, 1, 1, 1, 1, 1, 0,
  0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 0, 1, 1, 1, 0, 1, 0,

  1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1, 1, 0,
  0, 0, 1, 1, 0, 1, 1, 1, 1, 0,
  0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
  0, 0, 1, 1, 0, 0, 1, 1, 0, 1,

  1, 1, 0, 0, 1, 1, 0, 0, 0, 0,
  1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
  1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
  1, 1, 0, 0, 1, 1, 0, 0, 0, 1,
  1, 1, 1, 0, 1, 1, 0, 0, 1, 1,
  1, 1, 0, 0, 1, 1, 0, 0, 1, 0,

  1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0
)

data <- matrix(as.integer(data),n,M,byrow=TRUE)

true_components <-
  c(1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5)


# FIND THE CONDITIONAL DISTRIBUTION FOR THE MIXTURE COMPONENT OF AN OBSERVATION.

cond <- function (s, i, cntx, cnt1)
{ p <- cntx+1
  for (x in 1:X)
  { l <- 1
    for (j in 1:M)
    { p1 <- (cnt1[x,j]+1) / (cntx[x]+2)
      if (data[i,j]==1)
      { l <- l * p1
      }
      else
      { l <- l * (1-p1)
      }
    }
    p[x] <- p[x]*l
  }
  p/sum(p)
}


# MIXTURE SAMPLING RUN.  Does K iterations, each of which consists of
# updates of sites starting with initial state 's', using function
# 'indexes' (taking number of variables as argument) to produce the
# vector of variable indexes to use in an iteration (must always be
# the same length, default sequential), and using the function
# 'method' to modify the Gibbs sampling probabilities (default, no
# modification).  The value returned is a list with vectors 'eq10',
# and 'eq30' being the numbers of observations with the same mixture
# component as observations 10 or 30, after each update, 'cl1' being
# the cluster id for observation 1, after each update, 'self' being
# the frequency of self transitions, and 'final' being the final
# state.  Prints the frequency of self transitions, and means of
# 'cl1', 'eq10', and 'eq30'.
#
# Note that the counts needed to compute the posterior probability
# are computed at the beginning, but then updated incrementally after
# each change.

mix_run <- function (s, K, method=function(p,v)p, indexes=function(n)1:n)
{ n <- length(s)
  ix <- indexes(n)
  cl1 <- integer(K*n)
  eq10 <- integer(K*n)
  eq30 <- integer(K*n)
  cntx <- integer(X)
  cnt1 <- matrix(0L,X,M)
  for (x in 1:X) 
  { for (i in 1:n)
    { if (s[i]==x)
      { cntx[x] <- cntx[x] + 1L
        for (j in 1:M)
        { cnt1[x,j] <- cnt1[x,j] + data[i,j]
        }
      }
    }
  }
  self_trans <- 0
  self_trans_pr <- 0
  min_self_trans_pr <- 0
  prhalf <- 0
  for (t in 1:K)
  { for (k in 1:n)
    { w <- ix[k]
      x <- s[w]
      cntx[x] <- cntx[x] - 1L
      for (j in 1:M) cnt1[x,j] <- cnt1[x,j] - data[w,j]
      pi <- cond(s,w,cntx,cnt1)
      mxpi <- max(pi)
      if (mxpi>=0.5) prhalf <- prhalf+1
      msp <- if (mxpi<=0.5) 0 else 2*mxpi-1
      min_self_trans_pr <- min_self_trans_pr + msp
      p <- method(pi,x)
      self_trans_pr <- self_trans_pr + p[x]
      new <- sample(length(p),1,prob=p)
      if (new==x)
      { self_trans <- self_trans + 1
      }
      else
      { s[w] <- new
      }
      x <- s[w]
      cntx[x] <- cntx[x] + 1L
      for (j in 1:M) cnt1[x,j] <- cnt1[x,j] + data[w,j]
      h <- (t-1)*n + k
      cl1[h] <- as.integer(s[1]==1)
      eq10[h] <- cntx[s[10]]
      eq30[h] <- cntx[s[30]]
    }
    ix <- indexes(n)
  }

  if (TRUE)
  { cat("cntx:\n"); print(cntx)
    cat("cnt1:\n"); print(cnt1)
  }

  if (typeof(cl1)!="integer" || typeof(eq10)!="integer" 
                             || typeof(eq30)!="integer")
  { cat("WARNING: function values not of type integer\n")
  }

  cat("self transition frequency:",self_trans/(K*n),"\n")
  cat("average self transition probability:",self_trans_pr/(K*n),"\n")
  cat("minimum self transition probability:",min_self_trans_pr/(K*n),"\n")
  cat("fraction of time max pr > 1/2:",prhalf/(K*n),"\n")

  cat("\ncl1: mean",mean(cl1)," variance",var(cl1),"\n")
  cat("\neq10: mean",mean(eq10)," variance",var(eq10),"\n")
  cat("\neq30: mean",mean(eq30)," variance",var(eq30),"\n\n")

  names(s) <- c("A","B","C","D","E")[true_components]
  list (cl1=cl1, eq10=eq10, eq30=eq30, 
        self=self_trans/(K*n), prhalf=prhalf/(K*n),
        self_pr=self_trans_pr/(K*n), min_self_pr=min_self_trans_pr/(K*n),
        final=s)
}

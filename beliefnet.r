# BELIEF NETWORK DISTRIBUTION.


# NETWORK ARCHITECTURE.

N0 <- 3    # Number of low-level nodes
N1 <- 5    # Number of middle-level nodes
N2 <- 2    # Number of high-level nodes

M0 <- 3    # Number of possible values for low-level nodes
M1 <- 4    # Number of possible values for middle-level nodes
M2 <- 5    # Number of possible values for high-level nodes

total_values <- M0^N0 * M1^N1 * M2^N2

n <- N0 + N1 + N2


# RANDOMLY GENERATE PARAMETERS FOR THE NETWORK.

net_params <- function (s, t)
{
  p0 <- array (s*rt(N0*N1*M0*M1,t), c(N0,N1,M0,M1))
  p1 <- array (s*rt(N1*N2*M1*M2,t), c(N1,N2,M1,M2))
  p2 <- array (s*rt(N2*M2,t), c(N2,M2))

  list(p0=p0,p1=p1,p2=p2)
}

set.seed(2); par <- net_params (1, 4)    # Random parameters of network
cat("\nNetwork parameters:\n")
print(par)
cat("\n")


# FIND JOINT PROBABILITY OF VALUES IN NETWORK.

net_prob <- function (val, par)
{ 
  pr <- 1

  for (k in 1:N2)
  { e <- exp(par$p2[k,])
    pr <- pr * e[val[N0+N1+k]]/sum(e)
  }

  for (j in 1:N1)
  { w <- 0
    for (k in 1:N2)
    { w <- w + par$p1[j,k,,val[N0+N1+k]]
    }
    e <- exp(w)
    pr <- pr * e[val[N0+j]]/sum(e)
  }

  for (i in 1:N0)
  { w <- 0
    for (j in 1:N1)
    { w <- w + par$p0[i,j,,val[N0+j]]
    }
    e <- exp(w)
    pr <- pr * e[val[i]]/sum(e)
  }

  pr
}


# RETURN M'TH POSSIBLE SET OF NODE VALUES (USING LEXICOGRAPHICAL ORDER).

net_value <- function (m)
{ 
  m <- m-1  # change to starting at zero
  val <- numeric(N0+N1+N2)

  for (i in 1:N0)
  { val[i] <- m %% M0 + 1
    m <- m %/% M0
  }
  
  for (j in 1:N1)
  { val[N0+j] <- m %% M1 + 1
    m <- m %/% M1
  }
  
  for (k in 1:N2)
  { val[N0+N1+k] <- m %% M2 + 1
    m <- m %/% M2
  }
  
  if (m!=0)
  { stop("m too big")
  }

  val
}


# FIND MARGINAL DISTRIBUTION FOR A NODE.

marginal0 <- function (par,i)
{ m <- numeric(M0)
  for (h in 1:total_values)
  { val <- net_value(h)
    pr <- net_prob(val,par)
    m[val[i]] <- m[val[i]] + pr
  }
  m
}

marginal1 <- function (par,j)
{ m <- numeric(M1)
  for (h in 1:total_values)
  { val <- net_value(h)
    pr <- net_prob(val,par)
    m[val[N0+j]] <- m[val[N0+j]] + pr
  }
  m
}

marginal2 <- function (par,k)
{ m <- numeric(M2)
  for (h in 1:total_values)
  { val <- net_value(h)
    pr <- net_prob(val,par)
    m[val[N0+N1+k]] <- m[val[N0+N1+k]] + pr
  }
  m
}


# FIND JOINT DISTRIBUTION FOR TWO NODES.

joint02 <- function (par,i,k)
{ J <- matrix(0,M0,M2)
  for (h in 1:total_values)
  { val <- net_value(h)
    pr <- net_prob(val,par)
    J[val[i],val[N0+N1+k]] <- J[val[i],val[N0+N1+k]] + pr
  }
  J
}


# FIND CONDITIONAL DISTRIUBTION FOR A NODE GIVEN OTHER NODE VALUES.

cond0 <- function (par, val, i)
{ 
  pr <- numeric(M0)

  for (v in 1:M0)
  { val[i] <- v
    pr[v] <- net_prob(val,par)
  }

  pr / sum(pr)
}

cond1 <- function (par, val, i)
{ 
  pr <- numeric(M1)

  for (v in 1:M1)
  { val[N0+i] <- v
    pr[v] <- net_prob(val,par)
  }

  pr / sum(pr)
}

cond2 <- function (par, val, i)
{ 
  pr <- numeric(M2)

  for (v in 1:M2)
  { val[N0+N1+i] <- v
    pr[v] <- net_prob(val,par)
  }

  pr / sum(pr)
}


# BELIEF NET SAMPLING RUN.  Does K iterations, each of which consists
# of updates of values using parameters 'par', starting with initial
# state 'val', using function 'indexes' (taking number of variables as
# argument) to produce the vector of variable indexes to use in an
# iteration (must always be the same length, default sequential), and
# using the function 'method' to modify the Gibbs sampling
# probabilities (default, no modification).

beliefnet_run <- function (par, val, K, method=function(p,v) p, 
                                        indexes=function(n) 1:n)
{ 
  n <- length(val)
  ix <- indexes(n)

  mar1 <- integer(K*n)
  mar2 <- integer(K*n)
  and02 <- integer(K*n)
  self_trans <- 0
  self_trans_pr <- 0
  min_self_trans_pr <- 0
  prhalf<- 0
  for (t in 1:K)
  { # cat("Start",t,":",ix,":",val,"\n")
    for (k in 1:n)
    { w <- ix[k]
      x <- val[w]
      # cat("Gibbs",k,w,x,":",val,"\n")
      pi <- ( if (w<=N0) cond0(par,val,w)
              else if (w<=N0+N1) cond1(par,val,w-N0)
              else cond2(par,val,w-N0-N1) )
      mxpi <- max(pi)
      if (mxpi>=0.5) prhalf <- prhalf+1
      msp <- if (mxpi<=0.5) 0 else 2*mxpi-1
      min_self_trans_pr <- min_self_trans_pr + msp
      p <- method(pi,x)
      self_trans_pr <- self_trans_pr + p[x]
      # cat("pi:",pi," p:",p,"\n")
      new <- sample(length(p),1,prob=p)
      if (new==x)
      { self_trans <- self_trans + 1
      }
      else
      { val[w] <- new
        # cat("Change",w,"from",x,"to",new,"with",pr,"/",p,"\n")
      }
      x <- val[w]
      h <- (t-1)*n + k
      mar1[h] <- as.integer (val[N0+1]==1)
      mar2[h] <- as.integer (val[N0+N1+1]==1)
      and02[h] <- as.integer (val[1]==1 & val[N0+N1+1]==1)
    }
    ix <- indexes(n)
    # cat("End",t,":",val,"\n")
  }

  if (typeof(mar1)!="integer" || typeof(and02)!="integer")
  { cat("WARNING: function values not of type integer\n")
  }

  cat("self transition frequency:",self_trans/(K*n),"\n")
  cat("average self transition probability:",self_trans_pr/(K*n),"\n")
  cat("minimum self transition probability:",min_self_trans_pr/(K*n),"\n")
  cat("fraction of time max pr >= 1/2:",prhalf/(K*n),"\n")

  cat("\nmar1: mean",mean(mar1)," variance",var(mar1),"\n")
  cat("\nmar2: mean",mean(mar2)," variance",var(mar2),"\n")
  cat("\nand02: mean",mean(and02)," variance",var(and02),"\n")

  list (mar1=mar1, mar2=mar2, and02=and02,
        self=self_trans/(K*n), prhalf=prhalf/(K*n),
        self_pr=self_trans_pr/(K*n), min_self_pr=min_self_trans_pr/(K*n),
        final=val)
}

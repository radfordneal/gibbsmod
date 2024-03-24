# METHODS FOR MODIFYING GIBBS SAMPLING TRANSITION PROBABILITIES.
#
# All trans_... methods are passed the Gibbs sampling probabilities 
# (ie, conditional probabilities given other variables), pi, and the 
# current value, k, and return the vector of transition probabilities 
# to the various values.
#
# All sample_... methods are passed Gibbs sampling probabilities, the
# current value, and two uniform random variates, u1 and u2, in [0,1], 
# and return a new value sampled according to the modified transition
# probabilities using these random variates (possibly only one of them).
#
# Some of these methods are also passed an ordering of values, o.


# SAMPLE A VALUE ACCORDING TO A PROBABILITY VECTOR AND UNIFORM VARIATE.
#
# Will never return a value that has zero probability.

sample_value <- function (p,u)
{ 
  s <- 0
  for (i in 1:length(p))
  { if (p[i]>0)
    { s <- s+p[i]
      j <- i
      if (u<s)
      { break
      }
    }
  }

  return (j)
}


# ESTIMATE TRANSITION PROBABILITIES FROM SAMPLING, WITH N STRATIFIED CALLS.

trans_from_sample <- function (N, sample, pi, k, ...)
{
  p <- numeric(length(pi))

  for (i in 1:N)
  { u1 <- (i-1+runif(1)) / N
    u2 <- runif(1)
    j <- sample(pi,k,u1,u2,...)
    p[j] <- p[j]+1
  }

  p/N
}


# STANDARD GIBBS SAMPLING (GS).  Returns the conditional probabilities
# unchanged, ignoring the current value.

trans_GS <- function (pi,k) pi


# SAMPLE WITH GS.

sample_GS <- function (pi,k,u1,u2) sample_value (pi,u1)


# METROPOLIS HASTINGS GIBBS SAMPLING (MHGS).

trans_MHGS.1 <- function (pi,k)   # Algorithm 1 in the paper.
{
  if (any (1-pi <= 0))   # avoid division by zero below, reverting to Gibbs
  { return(pi)           #   sampling (with one probability almost 1)
  }

  p <- pmin (1, pi/(1-pi[k]), pi/(1-pi))  # min with 1 guards against rounding.

  p[k] <- 0
  p[k] <- max (0, 1-sum(p))  # avoid possible negative prob from round-off error

  p
}

trans_MHGS.2 <- function (pi,k)
{
  # Find proposal probabilities, v.  If all values other than k have zero
  # probability, they will be undefined, and we just return pi.

  v <- pi
  v[k] <- 0
  s <- sum(v)
  if (s==0)
  { return(pi)
  }
  v <- v / s

  # Reduce proposal probabilities by probabability of acceptance, while
  # accumulating self transition probability from rejection probabilities.

  for (j in 1:length(pi))
  { if (j!=k && pi[j]<pi[k])
    { a <- (1-pi[k])/(1-pi[j])
      v[k] <- v[k]+v[j]*(1-a)
      v[j] <- v[j]*a
    }
  }

  v
}

trans_MHGS <- trans_MHGS.1


# SAMPLE WITH MHGS.

sample_MHGS <- function (pi,k,u1,u2)
{
  p <- pi
  p[k] <- 0
  s <- sum(p)

  if (s==0)
  { return(k)
  }

  i <- sample_value (p/s, u1)

  if (u2 < (1-pi[k])/(1-pi[i])) i else k
}


# NESTED ANTITHETIC MODIFICATION METHOD (NAM).

trans_NAM.1 <- function (pi,k,o=1:m)   # Algorithm 2 in the paper
{
  m <- length(pi)

  epsilon <- 0   # could instead be set to a small non-zero value

  p <- rep(NA,m)

  s <- 1  # sum of probabilities for values that have not yet been focal
  f <- 1  # sum of transition prob from k to values that have not yet been focal

  # Find modified transition probabilities from the current value to 
  # successive focal values, until the focal value is the current value.

  i <- 1

  while (o[i] != k)
  { 
    # After seeing a focal value with probability at least as large as
    # remaining values, just store zeros.  (If epsilon to greater than 0
    # above, avoids tiny probabilities from rounding.)

    if (f <= epsilon)
    { p[o[i]] <- 0
    }

    else
    {
      # Let q be the original probability of the focal value; update s to be
      # the sum of probabilities for remaining non-focal values.

      q <- pi[o[i]]
      s <- s-q

      # Compute the transition probability from current value to focal value,
      # and find the new total probability for transitions to remaining values.

      if (q >= s)
      { p[o[i]] <- f
        f <- 0
      }
      else
      { p[o[i]] <- (q / s) * f  # will have p[o[i]]<=f<=1, even with rounding
        f <- f - p[o[i]]
      }
    }

    # cat("i:",i,"s:",round(s,7),"\n")
    # cat("i:",i+1,"f:",round(f,7),"\n")

    i <- i + 1
  }

  # Compute modified transition probabilities from the current value, k, which
  # is now focal, to values that have not previously been focal, as well as 
  # the self transition probability for k.

  if (f <= epsilon)
  { p[k] <- 0
    if (i<m) p[o[(i+1):m]] <- 0
  }
  else
  { q <- pi[k]
    s <- s-q
    if (q > s)
    { p[k] <- ((q-s) / q) * f  # will have p[k]<=f<=1, even with rounding
      if (i<m) p[o[(i+1):m]] <- pmin (f, (pi[o[(i+1):m]] / q) * f)
    }
    else
    { p[k] <- 0
      if (i<m) p[o[(i+1):m]] <- pmin (f, (pi[o[(i+1):m]] / s) * f)
    }
  }

  p
}

trans_NAM.2 <- function (pi,k,o=1:m)   # Algorithm 2 in the paper, tweaked for R
{
  m <- length(pi)

  epsilon <- 0      # could instead be set to a small non-zero value

  p <- numeric(m)   # transition probabilities are initially all zero

  s <- 1  # sum of probabilities for values that have not yet been focal
  f <- 1  # sum of transition prob from k to values that have not yet been focal

  # Find modified transition probabilities from the current value to 
  # successive focal values, until the focal value is the current value.

  i <- 1

  while (o[i] != k)
  { 
    # Let q be the original probability of the focal value; update s to be
    # the sum of probabilities for remaining non-focal values.

    q <- pi[o[i]]
    s <- s-q

    # Compute the transition probability from current value to focal value,
    # and find the new total probability for transitions to remaining values.

    if (q >= s)
    { p[o[i]] <- f
      f <- 0
    }
    else
    { p[o[i]] <- (q / s) * f  # will have p[o[i]]<=f<=1, even with rounding
      f <- f - p[o[i]]
    }

    # When f is zero/small, we can just leave the rest of p at zero and return.

    if (f <= epsilon)
    { return(p)
    }

    i <- i + 1
  }

  # Compute modified transition probabilities from the current value, k, which
  # is now focal, to values that have not previously been focal, as well as
  # the self transition probability for k.

  q <- pi[k]
  s <- s-q
  if (q > s)
  { p[k] <- ((q-s) / q) * f  # will have p[k]<=f<=1, even with rounding
    if (i<m) p[o[(i+1):m]] <- pmin (f, (pi[o[(i+1):m]] / q) * f)
  }
  else
  { p[k] <- 0
    if (i<m) p[o[(i+1):m]] <- pmin (f, (pi[o[(i+1):m]] / s) * f)
  }

  p
}

trans_NAM <- trans_NAM.2


# UPWARD NESTED ANTITHETIC MODIFICATION (UNAM).

trans_UNAM.0 <- function (pi,k) trans_NAM(pi,k,order(pi))

trans_UNAM.1 <- function (pi,k)  # Algorithm 3 in the paper
{
  m <- length(pi)
  o <- order(pi)

  p <- rep(NA,m)

  s <- 1  # sum of probabilities for values that have not yet been focal
  f <- 1  # sum of transition prob from k to values that have not yet been focal

  # Find modified transition probabilities from the current value to 
  # successive focal values, until the focal value is the current value.

  i <- 1

  while (o[i] != k)
  { 
    # Let q be the original probability of the focal value; update s to be
    # the sum of probabilities for remaining non-focal values.

    q <- pi[o[i]]
    s <- s-q

    # Compute the transition probability from current value to focal value,
    # and find the new total probability for transitions to remaining values.

    p[o[i]] <- pmin (f, (q / s) * f)  # min with f in case q>s due to rounding
    f <- f - p[o[i]]

    i <- i + 1
  }

  # Compute modified transition probabilities from the current value, k, which
  # is now focal, to values that have not previously been focal, as well as
  # the self transition probability for k.

  if (i==m)
  { p[k] <- f
  }
  else
  { p[k] <- 0
    p[o[(i+1):m]] <- pmin (f, (pi[o[(i+1):m]] / (s-pi[k])) * f)
  }

  p
}

trans_UNAM.iter <- function (pi,k,times=100)  # iterate MH modification
{
  m <- length(pi)
  o <- order(pi)

  T <- matrix(pi,m,m,byrow=TRUE)

  for (i in 1:times)
  { D <- diag(T)
    if (any(1-D<=0)) break
    T <- pmin (T / matrix(1-D,m,m,byrow=TRUE),
               T / matrix(1-D,m,m,byrow=FALSE))
    diag(T) <- 0
    diag(T) <- 1-rowSums(T)
  }

  T[k,]
}

trans_UNAM <- trans_UNAM.1     # Which algorithm to use


# DOWNWARD NESTED ANTITHETIC MODIFICATION (DNAM).

trans_DNAM.0 <- function (pi,k) trans_NAM(pi,k,o=rev(order(pi)))

trans_DNAM.1 <- function (pi,k)   # Algorithm 4 in the paper
{
  # Quickly handle the case where current value has probability half or more,
  # without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  # Otherwise call the general trans_NAM procedure, giving reverse order by
  # probability.  Don't use order with decreasing=TRUE, since for consistency
  # in comparisons it's desirable for the order be the reverse of UNAM when
  # there are ties in probability.

  trans_NAM(pi,k,o=rev(order(pi)))
}

trans_DNAM <- trans_DNAM.1


# AVERAGE OF UNAM AND DNAM.

trans_UDNAM <- function (pi,k) (trans_UNAM(pi,k) + trans_DNAM(pi,k)) / 2


# FIND THE DISTRIBUTION FOR THE NEXT STATE WITH ZERO-SELF DNAM.

trans_ZDNAM.1 <- function (pi,k)    # Algorithm 5 in the paper
{

  # Quickly handle the case where current value has probability half or more,
  # without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  m <- length(pi)
  p <- rep(NA,m)

  o <- rev(order(pi))  # reverse prob order; not done with decreasing=TRUE
                       # to simplify comparison with other methods when ties

  # Handle case where a value has probability of 1/2 or more.  Won't be
  # the current value, since that's handled above.

  if (pi[o[1]] >= 0.5)
  { p[o[1]] <- 1
    for (i in 2:m)
    { p[o[i]] <- 0
    }
    return(p)
  }

  s <- 1  # sum of probabilities for values that have not yet been focal
  f <- 1  # sum of transition prob from k to values that have not yet been focal

  # Find modified transition probabilities from the current value to
  # successive focal values, until the focal value is the current value, or
  # special handling to avoid a non-zero self-transition probability is needed.

  i <- 1

  while (f > 0 && o[i] != k && pi[o[i+1]] < s-pi[o[i]]-pi[o[i+1]])
  {
    # Let q be the probability of the focal value; update s to be
    # the sum of probabilities for remaining non-focal values.

    q <- pi[o[i]]
    s <- s-q                # guaranteed to be positive, even with rounding

    # Compute the transition probability from current value to focal value,
    # and find the new total probability for transitions to remaining values.

    p[o[i]] <- f * (q / s)  # guaranteed p[o[i]] <= f <= 1, even with rounding
    f <- f - p[o[i]]

    i <- i + 1
  }

  q <- pi[o[i]]
  s <- s-q

  if (f > 0 && s > 0 && i < m)
  {
    q2 <- pi[o[i+1]]
    s2 <- max (0, s-q2)  # max guards against round-off error
  
    if (q2 >= s2)
    { 
      # Use the special construction to avoid a non-zero self-transition 
      # probability.
  
      A <- (q+q2-s2) / 2
      if (k == o[i])
      { p[o[i]] <- 0
        p[o[i+1]] <- f * A / q
      }
      else if (k == o[i+1])
      { p[o[i]] <- f * A / q2
        p[o[i+1]] <- 0
      }
  
      if (s2<=0)
      { i <- i+2
      }
      else
      { B <- (q-q2+s2) / (2*s2)
        C <- (s2+q2-q) / (2*s2)
        if (k == o[i])
        { i <- i+2
          while (i <= m)
          { p[o[i]] <- f * B * pi[o[(i)]] / q
            i <- i+1
          }
        }
        else if (k == o[i+1])
        { i <- i+2
          while (i <= m)
          { p[o[i]] <- f * C * pi[o[i]] / q2
            i <- i+1
          }
        }
        else
        { p[o[i]] <- f * B
          p[o[i+1]] <- f * C
          i <- i+2
        }
      }
    }
    else
    {
      # Compute modified transition probabilities from the current value 
      # (now focal) to values that have not previously been focal.

      p[o[i]] <- 0
      i <- i+1
      while (i <= m)
      { p[o[i]] <- (pi[o[i]] / s) * f
        i <- i+1
      }
    }
  }

  # Set any remaining transition probabilities to zero.

  while (i <= m)
  { p[o[i]] <- 0
    i <- i+1
  }

  return(p)
}

trans_ZDNAM.2 <- function (pi,k)    # Algorithm 5 in the paper, tweaked for R
{
  # Quickly handle the case where current value has probability half or more,
  # without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  m <- length(pi)

  o <- rev(order(pi))  # reverse prob order; not done with decreasing=TRUE
                       # to simplify comparison with other methods when ties

  p <- numeric(m)      # transition probabilities are initially all zero

  # Handle case where a value has probability of 1/2 or more.  Won't be
  # the current value, since that's handled above.

  if (pi[o[1]] >= 0.5)
  { p[o[1]] <- 1
    return(p)
  }

  epsilon <- 0           # could instead be set to a small non-zero value

  s <- 1  # sum of probabilities for values that have not yet been focal
  f <- 1  # sum of transition prob from k to values that have not yet been focal

  # Find modified transition probabilities from the current value to
  # successive focal values, until the focal value is the current value, or
  # special handling to avoid a non-zero self-transition probability is needed.

  i <- 1

  while (o[i] != k && pi[o[i+1]] < s-pi[o[i]]-pi[o[i+1]])
  {
    # Let q be the probability of the focal value; update s to be
    # the sum of probabilities for remaining non-focal values.

    q <- pi[o[i]]
    s <- s-q                # guaranteed to be positive, even with rounding

    # Compute the transition probability from current value to focal value,
    # and find the new total probability for transitions to remaining values.

    p[o[i]] <- f * (q / s)  # guaranteed <= f <= 1, even with rounding
    f <- f - p[o[i]]

    # When f is zero/small, we can just leave the rest of p at zero and return.

    if (f <= epsilon)
    { return(p)
    }

    i <- i + 1
  }

  if (i < m)
  {
    q <- pi[o[i]]
    q2 <- pi[o[i+1]]

    s <- s-q
    s2 <- max (0, s-q2)  # max guards against round-off error
  
    if (q2 >= s2)
    { 
      # Use the special construction to avoid a non-zero self-transition 
      # probability.
  
      A <- (q+q2-s2) / 2
      if (k == o[i])
      { p[o[i+1]] <- f * A / q
      }
      else if (k == o[i+1])
      { p[o[i]] <- f * A / q2
      }
  
      if (s2>0)    
      { B <- (q-q2+s2) / (2*s2)
        C <- (s2+q2-q) / (2*s2)
        if (k == o[i])
        { if (i+1<m) p[o[(i+2):m]] <- f * B * pi[o[(i+2):m]] / q
        }
        else if (k == o[i+1])
        { if (i+1<m) p[o[(i+2):m]] <- f * C * pi[o[(i+2):m]] / q2
        }
        else
        { p[o[i]] <- f * B
          p[o[i+1]] <- f * C
        }
      }
    }
    else if (s > 0)
    {
      # Compute modified transition probabilities from the current value 
      # (now focal) to values that have not previously been focal.
  
      p[o[(i+1):m]] <- (pi[o[(i+1):m]] / s) * f
    }
  }

  return(p)
}

trans_ZDNAM <- trans_ZDNAM.2


# FIND THE DISTRIBUTION FOR THE NEXT STATE VARIABLE USING THE ST METHOD.
#
# Takes the shift amount (default, max pi) and and an ordering (default,
# original order) as optional arguments.
#
# Due to Suwa and Todo (2010).  Method from Suwa (2022), with correction.

trans_ST <- function (pi,k,s=max(pi),o=1:m)   # Algorithm 6 in the paper
{
  m <- length(pi)

  # Quickly handle the case where the current value has probability half or 
  # more, without needing to find cumulative probabilities.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  # Compute cumulative probabilities, in the order give by o, but stored
  # in the original order.  Set S to sum of all probabilities, which should
  # be one, but may differ due to rounding.

  C <- rep(NA,m)
  S <- 0
  for (i in 1:m)
  { C[o[i]] <- S
    S <- S + pi[o[i]]
  }

  # Find the flows from the current value to each value, and the total flow.

  delta <- pi[k] - s + C[k] - C   # delta[k] will be exactly zero if s=pi[k]
  delta2 <- delta + S

  v <- pmax (0, pmin (delta, pi[k] + pi - delta, pi[k], pi)) +
       pmax (0, pmin (delta2, pi[k] + pi - delta2, pi[k], pi))

  t <- sum(v)

  # If the total flow is zero, return transition probabilities that give 
  # probability 1 of moving to the most probable value.

  if (t==0)
  { p <- numeric(m)
    p[which.max(pi)] <- 1
    return(p)
  }

  # Return transition probabilities found by normalizing flows by their sum,
  # which should be pi(k), but may differ due to rounding.

  v/t
}


# SAMPLE WITH ST.  Written to use the original order only.

sample_ST.0 <- function (pi,k,u1,u2,s=max(pi))
{ 
  u <- if (k==1) 0 else sum(pi[1:(k-1)])
  u <- u + (u1*pi[k] - s)
  if (u<=0) u <- u+1
  sample_value(pi,u)
}


# SAMPLE WITH ST.  Algorithm 7 in paper.

sample_ST <- function (pi, k, u1, u2, o=1:m, s=max(pi))
{ 
  m <- length(pi)

  # Find the sum, u, of probabilities of values before k in the ordering used.

  i <- 1
  u <- 0
  while (o[i] != k)
  { u <- u + pi[o[i]]
    i <- i+1
  }

  # Add a random amount to u, while subtracting the shift, with wrap-around.

  u <- u + (u1*pi[k] - s)  # guaranteed not to increase u when s=max(pi)
  if (u<=0) u <- u+1       # guarantees that u is greater than zero

  # Use this value of u to pick a value to transition to, picking an arbitrary
  # value with non-zero probability if no value chosen due to round-off error.

  i <- 0
  s <- 0
  while (i < m && u > s)
  { i <- i+1
    if (pi[o[i]]>0)
    { s <- s+pi[o[i]]
      j <- o[i]
    }
  }

  return(j)
}


# FIND THE DISTRIBUTION FOR THE NEXT STATE VARIABLE USING DOWNWARD ST.

trans_DST <- function (pi,k)
{
  m <- length(pi)

  # Quickly handle the case where the current value has probability half or
  # more, without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  # Find order that puts values in non-increasing order of probability.
  # Highest-probability value will be at the bottom.

  o <- rev(order(pi)) # don't use decreasing=TRUE, so order same as in trans_UST
  pm <- pi[o[1]]

  # Specially handle case where a value has probability of half or more,
  # but not the current value, which would already have been handled above.

  if (pm>=0.5)
  { p <- numeric(m)
    p[o[1]] <- 1
    return(p)
  }

  # Return the result of using the general ST algorithm.

  trans_ST (pi, k, s=pm, o=o)
}


# FIND THE DISTRIBUTION FOR THE NEXT STATE VARIABLE USING UPWARD ST.

trans_UST <- function (pi,k)
{
  m <- length(pi)

  # Quickly handle the case where the current value has probability half or
  # more, without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  # Find order that puts values in non-decreasing order of probability.
  # Highest-probability value will be at the top.

  o <- order(pi)
  pm <- pi[o[m]]

  # Specially handle case where a value has probability of half or more,
  # but not the current value, which would already have been handled above.

  if (pm>=0.5)
  { p <- numeric(m)
    p[o[m]] <- 1
    return(p)
  }

  # Return the result of using the general ST algorithm.

  trans_ST (pi, k, s=pm, o=o)
}


# REVERSIBLE METHOD OBTAINED BY AVERAGING UST AND DST TRANSITIONS.

trans_UDST <- function (pi,k) (trans_UST(pi,k) + trans_DST(pi,k)) / 2


# FIND THE DISTRIBUTION FOR THE NEXT STATE VARIABLE USING HST.
#
# Due to Suwa (2022).

trans_HST.1 <- function (pi,k)
{
  trans_ST (pi, k, s=0.5)
}

trans_HST.2 <- function (pi,k)
{
  m <- length(pi)
  p <- numeric(m)

  # Handle transition from a zero-probability value specially, to avoid
  # division by zero below.  What is done does not affect correctness.

  if (pi[k]<=0)
  { p[which.max(pi)] <- 1
    return(p)
  }

  # Quickly handle the case where current state has probability half or more,
  # without needing to search for value with maximum probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  # Find value with maximimum probability.

  w <- which.max(p)
  pm <- pi[w]

  # Specially handle case where a state has probability of 1/2 or more
  # (not the current one, which would already have been handled above).

  if (pm>=0.5)
  { p[w] <- 1
    return(p)
  }

  # Set o to order that puts the max probability value at end (and doesn't
  # change order otherwise).

  o <- c(seq_len(w-1),w+seq_len(m-w),w)

  # Find the indexes of the current value and the value where the cumulative
  # probability first crosses 1/2.

  i <- m
  s <- 0
  i_k <- NULL
  i_half <- NULL
  while (is.null(i_k) || is.null(i_half))
  { s <- s + pi[o[i]]
    if (o[i] == k)
    { i_k <- i
      s_k <- s
    }
    if (is.null(i_half) && s >= 0.5)
    { i_half <- i
      s_half <- s
    }
    i <- i - 1
  }

  # Handle situation where the region for the current value extends above 1/2.

  if (s_k > 0.5)
  { j <- m
    s <- 0
    while (s+0.5 < s_k)
    { t <- s+pi[o[j]]
      z <- min(s_k,t+0.5) - max(s_k-pi[k],s+0.5)
      if (z > 0)
      { p[o[j]] <- z / pi[k]
      }
      s <- t
      j <- j - 1
    }
  }

  # Handle situation where the region for the current value starts below 1/2.

  if (s_k-pi[k] < 0.5)
  { j <- i_half
    s <- s_half - pi[o[i_half]]
    while (j > 0) # && s < s_k)
    { t <- s + pi[o[j]]
      z <- min(s_k,t-0.5) - max(s_k-pi[k],s-0.5)
      if (z > 0)
      { p[o[j]] <- z / pi[k]
      }
      s <- t
      j <- j - 1
    }
  }

  return(p)
}

trans_HST <- trans_HST.1


# FIND THE DISTRIBUTION FOR THE NEXT STATE VARIABLE USING ORDERED HST.

trans_OHST.1 <- function (pi,k)
{
  # Quickly handle the case where the current value has probability half or
  # more, without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  trans_ST (pi, k, s=0.5, o=order(pi))
}

trans_OHST.2 <- function (pi,k)
{
  m <- length(pi)
  p <- numeric(m)

  # Handle transition from a zero-probability value specially, to avoid
  # division by zero below.  What is done does not affect correctness.

  if (pi[k]<=0)
  { p[which.max(pi)] <- 1
    return(p)
  }

  # Quickly handle the case where current state has probability half or more,
  # without needing to order values by probability.

  if (pi[k]>=0.5)
  { p <- pmin (1, pi / pi[k])     # min guards against round-off error
    p[k] <- (2*pi[k]-1) / pi[k]
    return(p)
  }

  # Find order that puts values in non-decreasing order of probability.
  # Values will be stacked in [0,1] in order with highest-probability at bottom.

  o <- order(pi)

  # Specially handle case where a state has probability of 1/2 or more
  # (not the current one, which would already have been handled above).

  if (pi[o[m]]>=0.5)
  { p[o[m]] <- 1
    return(p)
  }

  # Find the indexes of the current value and the value where the cumulative
  # probability (going from most to least probable) first crosses 1/2.

  i <- m
  s <- 0
  i_k <- NULL
  i_half <- NULL
  while (is.null(i_k) || is.null(i_half))
  { s <- s + pi[o[i]]
    if (o[i] == k)
    { i_k <- i
      s_k <- s
    }
    if (is.null(i_half) && s >= 0.5)
    { i_half <- i
      s_half <- s
    }
    i <- i - 1
  }

  # Handle situation where the region for the current value extends above 1/2.

  if (s_k > 0.5)
  { j <- m
    s <- 0
    while (s+0.5 < s_k)
    { t <- s+pi[o[j]]
      z <- min(s_k,t+0.5) - max(s_k-pi[k],s+0.5)
      if (z > 0)
      { p[o[j]] <- z / pi[k]
      }
      s <- t
      j <- j - 1
    }
  }

  # Handle situation where the region for the current value starts below 1/2.

  if (s_k-pi[k] < 0.5)
  { j <- i_half
    s <- s_half - pi[o[i_half]]
    while (j > 0) # && s < s_k)
    { t <- s + pi[o[j]]
      z <- min(s_k,t-0.5) - max(s_k-pi[k],s-0.5)
      if (z > 0)
      { p[o[j]] <- z / pi[k]
      }
      s <- t
      j <- j - 1
    }
  }

  return(p)
}

trans_OHST <- trans_OHST.1


# FLATTENED SLICE SAMPLING.

trans_FSS.1 <- function (pi,k,zero_self=FALSE)  # Algorithm 8 in the paper.
{
  m <- length(pi)

  # Handle transition from zero-probability value specially.

  if (pi[k]<=0)
  { return(pi)
  }

  # Find the index of the most probable value, and its probability.

  ix1 <- which.max(pi)
  pi1 <- pi[ix1]

  # Handle the case where a value has probability of half or more specially.

  if (pi1>=0.5 || m<=2)           # m<=2 is a guard against round-off error
  { if (k==ix1)
    { p <- pi / pi1
      p[k] <- (2*pi1-1) / pi1
    }
    else
    { p <- numeric(m)
      p[ix1] <- 1
    }
    return(p)
  }

  # Find the probability of the second most probable value, and the difference
  # from that of the most probable value.

  pi2 <- 0
  for (i in 1:m)
  { if (i != ix1 && pi[i] > pi2)
    { pi2 <- pi[i]
    }
  }

  # Find the index of the value before the most probable value, or if 
  # zero_self is TRUE, the index of the first value before the most probable
  # value which will block movement beyond it from encountering a piece of the
  # most probable value. 

  ixb <- ix1
  while (TRUE)
  { ixb <- if (ixb==1) m else ixb-1
    po <- (0.5-pi1) + (0.5-pi[ixb]) # Computing this way reduces round-off error
    f <- (pi1-pi2)/po               # Guaranteed in [0,1] even with rounding
    if (!zero_self || pi[ixb] >= f*pi2)
    { break;
    }
  }

  # Find the part of the flow due to distributing the difference in probability
  # between most probable and second-most probable values among other values.
  # Here, f is the factor by which to multiply probabilities of values besides
  # ix1 and ixb to get the part of ix1 flowing there.

  v <- rep(NA,m)
  for (i in 1:m)
  { if (k==ix1 && i!=ix1 && i!=ixb)
    { v[i] <- f*pi[i]
    }
    else
    { v[i] <- 0
    }
  }

  # Find flow due to slice movement.

  l <- 0                           # lower end of probability region to move
  u <- if (k==ix1) pi2 else pi[k]  # upper end of probbability region to move

  i <- k

  while (l < u)
  { 
    # Move i backwards, going from ix1 to ixb, from ixb to before ix1, and
    # skipping ixb when otherwise going back.

    if (i==ix1) 
    { i <- ixb
    }
    else 
    { if (i==ixb)
      { i <- if (ix1==1) m else ix1-1
      }
      else
      { i <- if (i==1) m else i-1
      }
      if (i==ixb)
      { i <- if (ixb==1) m else ixb-1
      }
    }

    # Add to flow from slice movement of [l,u] region, and update l and u.

    if (l < pi[i])
    { if (i!=ix1 && i!=ixb)
      { t <- min (u, f*pi[i])
        if (l < t)
        { v[ix1] <- v[ix1] + t - l 
          l <- t
        }
      }
      t <- min (u, pi[i])
      v[i] <- v[i] + t - l
      l <- t
    }
  }

  # Return transition probabilities derived from flow.

  v / pi[k]
}

trans_ZFSS.1 <- function (pi,k) trans_FSS.1(pi,k,zero_self=TRUE)

trans_FSS.2 <- function (pi,k,zero_self=FALSE)  # Algorithm 8, tweaked for R
{
  m <- length(pi)

  # Handle transition from zero-probability value specially.

  if (pi[k]<=0)
  { return(pi)
  }

  # Find the index of the most probable value, and its probability.

  ix1 <- which.max(pi)
  pi1 <- pi[ix1]

  # Handle the case where a value has probability of half or more specially.

  if (pi1>=0.5 || m<=2)           # m<=2 is a guard against round-off error
  { if (k==ix1)
    { p <- pi / pi1
      p[k] <- (2*pi1-1) / pi1
    }
    else
    { p <- numeric(m)
      p[ix1] <- 1
    }
    return(p)
  }

  # Find the probability of the second most probable value, and the difference
  # from that of the most probable value.

  q <- pi
  q[ix1] <- 0
  pi2 <- max(q)
  d <- pi1 - pi2

  # Find the index of the value before the most probable value, or if 
  # zero_self is TRUE, the index of the first value before the most probable
  # value which will block movement beyond it from encountering a piece of the
  # most probable value. 

  ixb <- ix1
  while (TRUE)
  { ixb <- if (ixb==1) m else ixb-1
    po <- (0.5-pi1) + (0.5-pi[ixb]) # Computing this way reduces round-off error
    f <- d/po                       # Guaranteed in [0,1] even with rounding
    if (!zero_self || pi[ixb] >= f*pi2)
    { break;
    }
  }

  # Find the part of the flow due to distributing the difference in probability
  # between most probable and second-most probable values among other values.
  # Here, f is the factor by which to multiply probabilities of values besides
  # ix1 and ixb to get the part of ix1 flowing there.

  q[ixb] <- 0    # q now equals pi except that entries at ix1 and ixb are zero

  if (k==ix1)
  { v <- f*q
  }
  else
  { v <- numeric(m)
  }

  # Find flow due to slice movement.

  l <- 0                           # lower end of probability region to move
  u <- if (k==ix1) pi2 else pi[k]  # upper end of probbability region to move

  i <- k

  while (l < u)
  { 
    # Move i backwards, going from ix1 to ixb, from ixb to before ix1, and
    # skipping ixb when otherwise going back.

    if (i==ix1) 
    { i <- ixb
    }
    else 
    { if (i==ixb)
      { i <- if (ix1==1) m else ix1-1
      }
      else
      { i <- if (i==1) m else i-1
      }
      if (i==ixb)
      { i <- if (ixb==1) m else ixb-1
      }
    }

    # Add to flow from slice movement of [l,u] region, and update l and u.

    if (l < pi[i])
    { t <- min (u, f*q[i])
      if (l < t)
      { v[ix1] <- v[ix1] + t - l 
        l <- t
      }
      t <- min (u, pi[i])
      v[i] <- v[i] + t - l
      l <- t
    }
  }

  # Return transition probabilities derived from flow.

  v / pi[k]
}

trans_ZFSS.2 <- function (pi,k) trans_FSS.2(pi,k,zero_self=TRUE)

trans_FSS <- trans_FSS.2
trans_ZFSS <- trans_ZFSS.2


# SAMPLE WITH FSS/ZFSS

sample_FSS <- function (pi,k,u1,u2,zero_self=FALSE)  # Algorithm 9 in the paper.
{
  m <- length(pi)

  # Handle transition from zero-probability value specially.

  if (pi[k]<=0)
  { s <- 0
    i <- 0
    while (i < m && u1 >= s)
    { i <- i+1
      if (pi[i]>0)
      { s <- s + pi[i]
        j <- i
      }
    }
    return (j)
  }

  # Find the index of the most probable value, and its probability.

  ix1 <- which.max(pi)
  pi1 <- pi[ix1]

  # Handle the case where a value has probability of half or more specially.

  if (pi1>=0.5 || m<=2)           # m<=2 is a guard against round-off error
  { if (k!=ix1)
    { return(ix1)
    }
    else
    { s <- 0
      i <- 0
      while (i < m && u1 >= s)
      { i <- i + 1
        if (pi[i]>0)
        { if (i==k)
          { s <- s + (2*pi1-1) / pi1
          }
          else
          { s <- s + pi[i]/pi1
          }
          j <- i
        }
      }
      return (j)
    }
  }

  # Find the probability of the second most probable value, and the difference
  # from that of the most probable value.  Set j to the index of this value,
  # as a fallback if there's no other non-zero value other than ix1.

  pi2 <- 0
  for (i in 1:m)
  { if (i != ix1 && pi[i] > pi2)
    { pi2 <- pi[i]
      j <- i
    }
  }

  # Find the index of the value before the most probable value, or if 
  # zero_self is TRUE, the index of the first value before the most probable
  # value which will block movement beyond it from encountering a piece of the
  # most probable value. 

  ixb <- ix1
  while (TRUE)
  { ixb <- if (ixb==1) m else ixb-1
    po <- (0.5-pi1) + (0.5-pi[ixb]) # Computing this way reduces round-off error
    f <- (pi1-pi2)/po               # Guaranteed in [0,1] even with rounding
    if (!zero_self || pi[ixb] >= f*pi2)
    { break;
    }
  }

  # Multiply u1 by pi[k] to get r, which is uniform from 0 to the height of 
  # the bar for the current value.

  r <- u1*pi[k]

  # If the transition is from the most probable value, and r is in the region
  # that should be distributed to one of the values other than that and ixb,
  # then select such a value to be returned. 

  if (k==ix1 && r>=pi2)
  { r <- r - pi2
    s <- 0
    i <- 0
    while (i < m && r >= s)
    { i <- i + 1
      if (i!=ix1 && i!=ixb && pi[i]>0)
      { s <- s + f*pi[i]
        j <- i
      }
    }
    return (j)
  }

  # Return a value that is transitioned to due to slice movement.

  l <- 0                           # lower end of probability region to move
  u <- if (k==ix1) pi2 else pi[k]  # upper end of probbability region to move

  i <- k
  s <- 0

  while (l < u && r >= s)
  { 
    # Move i backwards, going from ix1 to ixb, from ixb to before ix1, and
    # skipping ixb when otherwise going back.

    if (i==ix1) 
    { i <- ixb
    }
    else 
    { if (i==ixb)
      { i <- if (ix1==1) m else ix1-1
      }
      else
      { i <- if (i==1) m else i-1
      }
      if (i==ixb)
      { i <- if (ixb==1) m else ixb-1
      }
    }

    # Look at slice movement from [l,u] region, and update l and u.

    if (l < pi[i])
    { if (i!=ix1 && i!=ixb)
      { t <- min (u, f*pi[i])
        if (l < t)
        { s <- s + t - l
          j <- ix1
          l <- t
        }
      }
      if (r >= s)
      { t <- min (u, pi[i])
        s <- s + t - l
        j <- i
        l <- t
      }
    }
  }

  return (j)
}

sample_ZFSS <- function (pi,k,u1,u2) sample_FSS(pi,k,u1,u2,zero_self=TRUE)


# ONE ANTITHETIC MODIFICATION OF TRANSITION MATRIX.

mat_AM <- function (pi,P,A,B,
                    delta = if (pi_A>pi_B) pi_B/pi_A else pi_A/pi_B)
{
  m <- length(pi)
  if (!all (A %in% (1:m)) || !all (B %in% (1:m)) || any(A %in% B))
  { stop("Invalid A or B")
  }
  if (nrow(P)!=m || ncol(P)!=m)
  { stop("Invalid P")
  }

  pi_A <- sum(pi[A])
  pi_B <- sum(pi[B])

  if (delta<=0)
  { stop("Invalid delta")
  }

  R <- P
  
  for (a in A)
  { R[a,A] <- R[a,A] - delta * pi[A] * pi_B / pi_A
    R[a,B] <- R[a,B] + delta * pi[B]
  }

  for (b in B)
  { R[b,B] <- R[b,B] - delta * pi[B] * pi_A / pi_B
    R[b,A] <- R[b,A] + delta * pi[A]
  }

  if (any(R < -1e-10))
  { cat("Negative probability\n")
  }

  R
}


# X

mat_X <-  function (pi,P,i,delta)
{
  m <- length(pi)
  if (P[i,i]<delta)
  { stop("Invalid delta")
  }

  R <- P
  R[i,i] <- R[i,i]-delta
  f <- (1-R[i,i]) / (1-P[i,i])
  R[i,-i] <- f * R[i,-i]
  R[-i,i] <- f * R[-i,i]

  for (j in 1:m)
  { if (j!=i)
    { R[j,-i] <- R[j,-i] * (1-R[j,i]) / sum(R[j,-i])
    }
  }

  R
}


# Y

mat_Y <-  function (pi,P,i,
           delta = (P[i,i]-sum(P[1,1:(i-1)])) / (1 + pi[i]/pi[i+1] 
                     * P[1,i+1]/sum(P[1,1:(i-1)])/sum(P[i+1,1:(i-1)])))
{
  m <- length(pi)
  if (i<=1 || i>=m)
  { stop("Invalid i")
  }
  if (P[i,i]<delta)
  { stop("Invalid delta")
  }

  R <- P

  R[i,i] <- R[i,i]-delta
  R[i,i+1] <- R[i,i+1]+delta
  R[i+1,i] <- R[i+1,i]+delta*pi[i]/pi[i+1]

  f <- delta*pi[i]/pi[i+1]/sum(P[i+1,1:(i-1)])

  R[i+1,1:(i-1)] <- R[i+1,1:(i-1)]*(1-f)
  R[1:(i-1),i+1] <- R[1:(i-1),i+1]*(1-f)

  g <- delta*P[1,i+1]*pi[i]/pi[i+1]/sum(P[1,1:(i-1)])/sum(P[i+1,1:(i-1)])

  R[1:(i-1),1:(i-1)] <- R[1:(i-1),1:(i-1)]*(1+g)

  R
}

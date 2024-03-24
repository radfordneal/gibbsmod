# VARIOUS SCAN ORDERS FOR GIBBS OR MODIFIED GIBBS UPDATES.
#
# These all take the number, n, of variables as their argument, and return a
# vector of that length with indexes from {1,...,n}.


# RANDOM INDEXES.  Returns n indexes chosen randomly, with repetition possible.

scan_random <- function (n) sample(n,n,replace=TRUE)


# SEQUENTIAL ORDER.

scan_sequential <- function (n) 1:n


# SEQUENTIAL ORDER AFTER SOME PERMUTATION.  The permuation is in the global
# variable shuffle_order, which must be set elsewhere, and is intended to
# be fixed to the same order for all runs.

scan_shuffled_sequential <- function (n) shuffle_order


# RANDOM ORDER.  A random permutation chosen again each time it is called.

scan_random_order <- function (n) sample(n)


# RANDOM ORDER REPEATED FOUR TIMES.  A random permutation, chosen again 
# after every four uses.

rand_perm <- NULL
rand_perm_uses <- 4

scan_random_order_x4 <- function (n)
{ if (rand_perm_uses==4)
  { rand_perm <<- sample(n)
    rand_perm_uses <<- 0
  }
  rand_perm_uses <<- rand_perm_uses + 1
  rand_perm
}


# RANDOM DIRECTION.  Either 1,...,n or n,...,1, randomly chosen each time.

scan_random_direction <- function (n) if (runif(1)<0.5) 1:n else n:1


# CHECKERBOARD.  Assumes variables are in a square array.  Does all of one
# colour of squares, then all other colour.

scan_checkerboard <- function (n)
{ if (n < 2) stop("checkerboard must be at least 2x2")
  m <- as.integer(sqrt(n))
  if (m^2 != n) stop("checkerboard must be square")
  b <- rep (c(rep(c(TRUE,FALSE),length=m),rep(c(FALSE,TRUE),length=m)),length=n)
  c((1:n)[b],(1:n)[!b])
}

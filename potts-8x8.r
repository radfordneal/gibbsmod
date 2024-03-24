# POTTS MODEL, 8x8 VERSION.

source("potts.r")

NC <- 8         # Number of columns in site array
NR <- 8         # Number of rows in site array
n <- NR*NC      # Total number of site variables
NV <- 4L        # Number of possible values for a site variable (1..NV)
b <- 0.85       # Strength of bonding between neighboring sites

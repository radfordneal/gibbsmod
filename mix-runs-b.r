# DO THE MIXTURE MODEL EXPERIMENTAL RUNS - SET B.

source("methods.r")
source("scans.r")
source("mix.r")
source("plot.r")

rtype <- "b"

meth <- list (GS=trans_GS,
              MHGS=trans_MHGS,
              UNAM=trans_UNAM,
              DNAM=trans_DNAM,
              UDNAM=trans_UDNAM,
              ZDNAM=trans_ZDNAM
        )

scan <- list (Random=scan_random,
              `Shuffled Sequential`=scan_shuffled_sequential,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("mix-runs-tail.r")

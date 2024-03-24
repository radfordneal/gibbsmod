# DO THE BELIEFNET EXPERIMENTAL RUNS - SET B.

source("methods.r")
source("scans.r")
source("beliefnet.r")
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
              Sequential=scan_sequential,
              `Shuffled Sequential`=scan_shuffled_sequential,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("beliefnet-runs-tail.r")

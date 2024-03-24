# DO THE BELIEFNET EXPERIMENTAL RUNS - SET B.

source("methods.r")
source("scans.r")
source("beliefnet.r")
source("plot.r")

rtype <- "c"

meth <- list (UNAM=trans_UNAM,
              ZDNAM=trans_ZDNAM,
              ST=trans_ST,
              UDST=trans_UDST,
              FSS=trans_FSS,
              ZFSS=trans_ZFSS
        )

scan <- list (Random=scan_random,
              Sequential=scan_sequential,
              `Shuffled Sequential`=scan_shuffled_sequential,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("beliefnet-runs-tail.r")

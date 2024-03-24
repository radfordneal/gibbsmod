# DO THE BELIEFNET EXPERIMENTAL RUNS - SET B.

source("methods.r")
source("scans.r")
source("beliefnet.r")
source("plot.r")

rtype <- "a"

meth <- list (ST=trans_ST, 
              DST=trans_DST,
              UST=trans_UST,
              UDST=trans_UDST,
              HST=trans_HST,
              OHST=trans_OHST
        )

scan <- list (Random=scan_random,
              Sequential=scan_sequential,
              `Shuffled Sequential`=scan_shuffled_sequential,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("beliefnet-runs-tail.r")

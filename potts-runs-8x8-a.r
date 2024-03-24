# DO THE 8x8 POTTS EXPERIMENTAL RUNS - SET A.

source("methods.r")
source("scans.r")
source("potts-8x8.r")
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
              Checkerboard=scan_checkerboard,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("potts-runs-8x8-tail.r")

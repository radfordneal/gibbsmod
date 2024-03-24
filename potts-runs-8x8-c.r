# DO THE 8x8 POTTS EXPERIMENTAL RUNS - SET C.

source("methods.r")
source("scans.r")
source("potts-8x8.r")
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
              Checkerboard=scan_checkerboard,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("potts-runs-8x8-tail.r")

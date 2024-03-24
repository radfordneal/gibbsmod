# DO THE 8x8 POTTS EXPERIMENTAL RUNS - SET B.

source("methods.r")
source("scans.r")
source("potts-8x8.r")
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
              Checkerboard=scan_checkerboard,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("potts-runs-8x8-tail.r")

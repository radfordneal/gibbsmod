# DO THE MIXTURE MODEL EXPERIMENTAL RUNS - SET A.

source("methods.r")
source("scans.r")
source("mix.r")
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
              `Shuffled Sequential`=scan_shuffled_sequential,
              `Random order`=scan_random_order,
              `Random order x4`=scan_random_order_x4
        )

source("mix-runs-tail.r")

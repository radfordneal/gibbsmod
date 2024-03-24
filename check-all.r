# RUN CHECKS ON ALL METHODS.

source("methods.r")
source("check.r")

random_checks (trans_GS, reversible=TRUE)

random_checks (trans_MHGS.1, reversible=TRUE, same_as=trans_MHGS.2)
random_checks (trans_MHGS.2, reversible=TRUE, same_as=trans_MHGS.1)

random_checks (trans_NAM.1, reversible=TRUE, same_as=trans_NAM.2)
random_checks (trans_NAM.2, reversible=TRUE, same_as=trans_NAM.1)

random_checks (trans_UNAM.0, reversible=TRUE, same_as=trans_UNAM.1)
random_checks (trans_UNAM.1, reversible=TRUE, same_as=trans_UNAM.0)
#random_checks (trans_UNAM.iter, reversible=TRUE, same_as=trans_UNAM.1)

random_checks (trans_DNAM.0, reversible=TRUE, same_as=trans_DNAM.1)
random_checks (trans_DNAM.1, reversible=TRUE, same_as=trans_DNAM.0)

random_checks (trans_UDNAM, reversible=TRUE)

random_checks (trans_ZDNAM.1, zero_self=TRUE, reversible=TRUE, 
                              same_as=trans_ZDNAM.2)
random_checks (trans_ZDNAM.2, zero_self=TRUE, reversible=TRUE, 
                              same_as=trans_ZDNAM.1)

random_checks (trans_ST, zero_self=TRUE)

random_checks (trans_DST, zero_self=TRUE, reverse_of=trans_UST)
random_checks (trans_UST, zero_self=TRUE, reverse_of=trans_DST)
random_checks (trans_UDST, zero_self=TRUE, reversible=TRUE)

random_checks (trans_HST, zero_self=TRUE, reversible=TRUE)
random_checks (trans_OHST, zero_self=TRUE, reversible=TRUE)

random_checks (trans_FSS.1, same_as=trans_FSS.2)
random_checks (trans_FSS.2, same_as=trans_FSS.1)

random_checks (trans_ZFSS.1,zero_self=TRUE,same_as=trans_ZFSS.2)
random_checks (trans_ZFSS.2,zero_self=TRUE,same_as=trans_ZFSS.1)


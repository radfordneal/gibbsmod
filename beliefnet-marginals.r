# COMPUTE MARGINAL PROBABILITIES FOR BELIEFNET MODEL.

source("beliefnet.r")

#marg0 <- marginal0(par,1)
print(marg1 <- marginal1(par,1))
print(marg2 <- marginal2(par,1))
print(joint <- joint02(par,1,1))


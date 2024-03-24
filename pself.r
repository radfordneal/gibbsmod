# PRODUCE PLOTS OF SELF-TRANSITION PROBABILITIES IN SIMPLE SITUATION.
#
# Scenario is maximum probability for a value is p, with other values
# having much smaller probability.

g <- seq(0.001,0.999,length=999)

mins <- function (p) pmax(0,(2*p-1))  # minimum self-transition probability
gs <- function (p) p^2                # self-transition probability of GS

pdf("runs/pself.pdf",width=9,height=5)
par(mfrow=c(1,2))

plot (g, mins(g), col="red", lwd=2, mgp=c(2,0.7,0),
      xlab="largest probability", 
      ylab="self-transition probability", 
      type="l", xaxs="i", xlim=c(0,1), yaxs="i", ylim=c(0,1))
lines(g, gs(g), col="blue", lwd=2)

plot (g, (1-mins(g))/(1-gs(g)), lwd=2, mgp=c(2,0.7,0),
     xlab="largest probability", 
     ylab="ratio of non-self-transition probabilities", 
     type="l", xaxs="i", xlim=c(0,1), yaxs="i", ylim=c(1,1.4))

dev.off()

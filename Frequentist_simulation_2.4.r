setwd("/Users/laura/Desktop")

source("auxiliary_functions.R")
source("final_function.R")

library(progress)
library(pcalg) 
library(MASS) 
library(parallel) 
library(abind) 
library(ggplot2) 
library(igraph)

sim2 <- function(n, p, truecovm, trueCPDAGamat, x, y, x.lab, y.lab) {
  
  # generate data according to trueDAG
  dat <- mvtnorm::rmvnorm(n=n, mean=rep(0, p), sigma= unname(truecovm))
  colnames(dat) <- colnames(trueCPDAGamat)
  
  # estimate CPDAG
  score = new("GaussL0penObsScore", dat)
  estCPDAG = ges(score)
  estCPDAGamat = 1*wgtMatrix(estCPDAG$essgraph)
  dimnames(estCPDAGamat) = list(1:p, 1:p)
  
  # estimate possible causal effects using different variants of IDA
  # true CPDAG
  
  trueLocal   = (ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(trueCPDAGamat))) 
  trueOptimal = (ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(trueCPDAGamat), 
                        method="optimal"))
  truePath   = (causal_effect_path(p = p, CPDAGamat = t(trueCPDAGamat),
                                      lab = 1:p, X = dat, x.pos = x, y.pos = y, 
                                      x.lab = x.lab, y.lab =y.lab))
  
  
  # estimated CPDAG
  estLocal    = (ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(estCPDAGamat)))
  estOptimal  = (ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(estCPDAGamat), 
                        method="optimal"))
  estPath     = (causal_effect_path(p = p, CPDAGamat = t(estCPDAGamat),
                                       lab = colnames(dat), X = dat, x.pos = x, y.pos = y, 
                                       x.lab = x.lab, y.lab =y.lab))
  #### Means ####
  # True
  tLocal_avg = mean(trueLocal)
  tOptimal_avg = mean(trueOptimal)
  tPath_avg  = mean(truePath)
  
  # Estimated
  eLocal_avg   = mean(estLocal)
  eOptimal_avg = mean(estOptimal)
  ePath_avg    = mean(estPath)
  
  # minimum absolute values
  # True
  tLocal_lb   = min(abs(trueLocal))
  tOptimal_lb = min(abs(trueOptimal))
  tPath_lb    = min(abs(truePath))
  
  eLocal_lb   = min(abs(estLocal))
  eOptimal_lb = min(abs(estOptimal))
  ePath_lb    = min(abs(estPath))
  
  return(cbind(tLocal_avg, tOptimal_avg, tPath_avg,
           eLocal_avg, eOptimal_avg, ePath_avg,
           tLocal_lb, tOptimal_lb, tPath_lb,
           eLocal_lb, eOptimal_lb, ePath_lb))
}



sim <- function(p, d, n){
  
  nonamenable <- FALSE
  nonnull <- FALSE
  while (!nonamenable | !nonnull) {
    
    # create a dag
    wfun <- function(m) {-1^rbinom(m, 1, 0.5) * runif(m, 0.1, 1)}
    trueDAG <- randDAG(p, d, method="er", wFUN=wfun)
    truecovm <- trueCov(trueDAG, back.compatible=TRUE)
    DAGamat <- t(as(trueDAG, "matrix"))
    
    # draw x and y
    ok <- FALSE
    while (!ok) {
      x <- sample(1:p, 1)
      x.lab = as.character(x)
      desX <- setdiff(possDe(DAGamat, x, possible=FALSE, type="dag"), x)
      if (length(desX) < 1) {next}
      if (length(desX)==1) {
        y <- desX
        y.lab = as.character(y)} 
      else {
        y <- sample(desX, 1)
        y.lab = as.character(y)}
      ok <- TRUE
    }
    
    # find the true CPDAG
    trueCPDAG <- dag2cpdag(trueDAG)
    trueCPDAGamat <- t(as(trueCPDAG, "matrix"))# NB: not encoded in our format! 
    dimnames(trueCPDAGamat) <- list(1:p, 1:p)
    
    nonamenable <- !pcalg:::isAmenable(trueCPDAGamat, x, y, type="cpdag")
    
    # find all true possible causal effects that can be estimated from the markov equiv class
    if (nonamenable) {
      # determine true smallest possible effect in CPDAG
      alltrue <- ida(x, y, trueCov(trueDAG, back.compatible=TRUE), t(trueCPDAGamat), 
                     type="pdag")
      true_LB  = min(abs(alltrue))
      true_avg = mean(alltrue)
      if (true_LB>(10^(-7))) {nonnull <- TRUE}
    }
  }
  
  # for 100 times simulate data and estimate the min abs of causal effect
  # using as input both the true cpdag (uncertainty comes just from data) and
  # the estimated cpdag (model uncertainty)
  res2 <- replicate(100, sim2(n , p, truecovm, trueCPDAGamat, x, y, x.lab, y.lab))
  
  # compute the MSE
  means = res2[ , 1:6, ] # est average causal effect according pa, opt, path in true cpdag and est cpdag
  lb    = res2[ , 7:12, ] # est lbs according to pa, opt, path, in true cpdag abd estimated cpdag
  
  
  MSE_avg_ce =  apply(means, 1, function (x) mean((x-true_avg)^2))
  MSE_lb_ce  =  apply(lb,    1, function (x) mean((x-true_LB)^2))
  
  res = c(unlist(MSE_avg_ce), unlist(MSE_lb_ce))
  names(res) <- c("MSE avg Pa True", "MSE avg Opt True", "MSE avg Path True", 
                  "MSE avg Pa Est", "MSE avg Opt Est", "MSE avg Path Est",
                  "MSE LB Pa True", "MSE LB Opt True", "MSE LB Path True",
                  "MSE LB Pa Est", "MSE LB Opt Est", "MSE LB Path Est")
  
  return(res)
  
  
}

nrep = 100

# Run the simulations:
# B -> p = 20, d =2, n = 1000
# A -> p = 100, d = 4, n = 1000



# Scenario B - n = 1000

pb <- progress_bar$new(format = "[:bar] :percent :elapsed Time Elapsed",
                       total = nrep)

set.seed(1234)
scenario_B_1000 <- replicate(nrep, {
  pb$tick()
  sim(p=20, d=2, n=1000)
})
pb$terminate() 

means = t(scenario_B_1000[1:6, ])
lbs   = t(scenario_B_1000[7:12, ])

write.csv(means, file = "means_B_freq_1000.csv", row.names = F)
write.csv(lbs, file = "lbs_B_freq_1000.csv", row.names = F)
# Scenario A

set.seed(1234)
pb <- progress_bar$new(format = "[:bar] :percent :elapsed Time Elapsed",
                       total = nrep)
scenario_A <- replicate(nrep, {
  pb$tick()
  sim(p=100, d=4, n=1000)
})
pb$terminate() 

means.A = t(scenario_A[1:6, ])
lbs.A   = t(scenario_A[7:12, ])

write.csv(means.A, file = "means_A_freq.csv", row.names = F)
write.csv(lbs.A, file = "lbs_A_freq.csv", row.names = F)


################
##### Plot #####
################


#1. Scenario A

#a. MSE of avg CE w.r.t. true causal effect 

means_A = read.csv("means_A_freq.csv")
lbs_A = read.csv("lbs_A_freq.csv")
geomean_A.means<- apply(means_A, 2, function(x) {exp(mean(log(x)))} )

# MSE avg CE w.r.t. true avg CE - true CPDAG (no model uncertainty)

true_A.means <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                                  each=100),
                            y=c(means_A[,1], means_A[,2], means_A[,3]) )

true_A.means$x <- factor(true_A.means$x, levels=c("IDA-parents", 
                                                  "IDA-Optimal", 
                                                  "f(L)-Path"))

# plot scenario A - true avg
ggplot(data=true_A.means, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_A.means[1]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_A.means[2]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_A.means[3]), shape=20, col ="red") + 
  ylim(0, 0.03) +
  annotate("text", label= "Scenario A - known CPDAG", x=2, y=0.025, size=4) +
  labs(x="", y="MSE of avg CE w.r.t. true avg CE") 


# plot scenario A - est avg

est_A.means <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                                 each=100),
                           y=c(means_A[,4], means_A[,5], means_A[,6]) )

est_A.means$x <- factor(est_A.means$x, levels=c("IDA-parents", 
                                                "IDA-Optimal", 
                                                "f(L)-Path"))

ggplot(data=est_A.means, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_A.means[4]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_A.means[5]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_A.means[6]), shape=20, col = "red") + 
  ylim(0, 0.05) +
  annotate("text", label=" Scenario A - Unknown CPDAG", x=2, y=0.05, size=4) +
  labs(x="", y="MSE of avg CE w.r.t. true avg CE") 

#b. Scenario A - MSE of LB w.r.t. true LB

geomean_A.lb<- apply(lbs_A, 2, function(x) {exp(mean(log(x)))} )


# MSE LB of CE - true CPDAG (no model uncertainty)

true_A.lb <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                               each=100),
                         y=c(lbs_A[,1], lbs_A[,2], lbs_A[,3]) )

true_A.lb$x <- factor(true_A.lb$x, levels=c("IDA-parents", "IDA-Optimal", "f(L)-Path"))

ggplot(data=true_A.lb, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_A.lb[1]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_A.lb[2]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_A.lb[3]), shape=20, col ="red") + 
  ylim(0, 0.02) +
  annotate("text", label="Scenario A - known CPDAG", x=2, y=0.02, size=4) +
  labs(x="", y="MSE of LB CE w.r.t. true LB CE") 


# est CPDAG lb

est_A.lb <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                              each=100),
                        y=c(lbs_A[,4], lbs_A[,5], lbs_A[,6]) )

est_A.lb$x <- factor(est_A.lb$x, levels=c("IDA-parents", "IDA-Optimal", "f(L)-Path"))

ggplot(data=est_A.lb, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_A.lb[4]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_A.lb[5]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_A.lb[6]), shape=20, col = "red") + 
  ylim(0, 0.05) +
  annotate("text", label=" Scenario A - unknown CPDAG", x=2, y=0.05, size=4) +
  labs(x="", y="MSE of estimated LB CE w.r.t. true LB CE")

# Table - Scenario A
apply(means_A, 2, mean)
apply(means_A, 2, sd)
apply(lbs_A,2, mean)
apply(lbs_A,2, sd)
apply(means_A, 2, median)
apply(lbs_A,2,median)
# Plot - Ratio of MSE - Scenario A

ratio_A = data.frame(x = rep(c("l:o MSE", "path:o MSE"), each = 100), 
                     y =c(est_A.lb[1:100,2] / est_A.lb[101:200,2],
                          est_A.lb[201:300,2] / est_A.lb[101:200,2]))

ratio_A$x <- factor(ratio_A$x, levels=c("l:o MSE", "path:o MSE"))

ggplot(data=ratio_A, aes(x=x, y=y)) +
  geom_violin(trim = F) + 
  geom_point(aes(x=1, y=exp(mean(log(ratio_A[1:100,2])))), shape=20, col = "red") + 
  geom_point(aes(x=2, y=exp(mean(log(ratio_A[101:200,2])))), shape=20, col = "red") + 
  geom_hline(yintercept=1, lty=2) +
  ylim(0,2) +
  annotate("text", label=" Scenario A - Estimated CPDAG",
           x=2, y=1.7, size=4) +
  labs(x="", y="RMSE")


#### Scenario B - n = 1000

means = read.csv("means_B_freq_1000.csv")
lbs = read.csv("lbs_B_freq_1000.csv")

geomean_B.means<- apply(means, 2, function(x) {exp(mean(log(x)))} )

# scenario b true avg CE
true_B.means <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                                  each=100),
                            y=c(means[,1], means[,2], means[,3]) )

true_B.means$x <- factor(true_B.means$x, levels=c("IDA-parents",
                                                  "IDA-Optimal",
                                                  "f(L)-Path"))

# plot scenario B - true avg
ggplot(data=true_B.means, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_B.means[1]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_B.means[2]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_B.means[3]), shape=20, col ="red") + 
  ylim(0, 0.008) +
  annotate("text", label=" Scenario B - known CPDAG", x=2, y=0.0075, size=4) +
  labs(x="", y="MSE of avg CE w.r.t. true avg CE") 

# plot scenario B - est avg

est_B.means <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                                 each=100),
                           y=c(means[ ,4], means[ ,5], means[ ,6]) )

est_B.means$x <- factor(est_B.means$x, levels=c("IDA-parents", "IDA-Optimal", 
                                                "f(L)-Path"))

ggplot(data=est_B.means, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_B.means[4]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_B.means[5]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_B.means[6]), shape=20, col = "red") + 
  ylim(0, 0.3) +
  annotate("text", label=" Scenario B - unknown CPDAG", x=2, y=0.3, size=4) +
  labs(x="", y="MSE of avg CE w.r.t. true avg CE") 

# the three methods behave the same 

# Scenario B - MSE of LB w.r.t. true LB

geomean_B.lb<- apply(lbs, 2, function(x) {exp(mean(log(x)))} )


# True LB
true_B.lb <- data.frame(x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                              each=100),
                        y=c(lbs[,1], lbs[,2], lbs[,3]) )

true_B.lb$x <- factor(true_B.lb$x, levels=c("IDA-parents", "IDA-Optimal", "f(L)-Path"))

# plot scenario B - true LB
ggplot(data=true_B.lb, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_B.lb[1]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_B.lb[2]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_B.lb[3]), shape=20, col ="red") + 
  ylim(0, 0.007) +
  annotate("text", label=" Scenario B - known CPDAG", x=2, y=0.006, size=4) +
  labs(x="", y="MSE of estimated LB CE w.r.t. true LB CE") 


# plot scenario B - est lb

est_B.lb <- data.frame( x=rep(c("IDA-parents", "IDA-Optimal", "f(L)-Path"),
                              each=100),
                        y=c(lbs[,4], lbs[,5], lbs[,6]) )

est_B.lb$x <- factor(est_B.lb$x, levels=c("IDA-parents", "IDA-Optimal", "f(L)-Path"))

ggplot(data=est_B.lb, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean_B.lb[4]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean_B.lb[5]), shape=20, col = "red") + 
  geom_point(aes(x=3, y=geomean_B.lb[6]), shape=20, col = "red") + 
  ylim(0, 0.40) +
  annotate("text", label=" Scenario B - unknown CPDAG", x=2, y=0.35, size=4) +
  labs(x="", y="MSE of estimated LB CE w.r.t. true LB CE")

# Table scenario B def

apply(means,2,mean)
apply(means,2,sd)
apply(means,2,median)
apply(lbs,2,mean)
apply(lbs,2,sd)
apply(lbs,2, median)

# Plot Ratio MSE - scenario B
ratio_B = data.frame(x = rep(c("l:o MSE", "path:o MSE"), each = 100), 
                     y =c(est_B.lb[1:100,2] / est_B.lb[101:200,2],
                          est_B.lb[201:300,2] / est_B.lb[101:200,2]))

ratio_B$x <- factor(ratio_B$x, levels=c("l:o MSE", "path:o MSE"))

ggplot(data=ratio_B, aes(x=x, y=y)) +
  geom_violin(trim = F) + 
  geom_point(aes(x=1, y=exp(mean(log(est_B.lb[1:100,2] / est_B.lb[101:200,2])))), 
             shape=20, col = "red") + 
  geom_point(aes(x=2, y=exp(mean(log(ratio_B[101:200,2])))), shape=20, col = "red") + 
  geom_point(aes(x=1, y= median(ratio_B[1:100,2])), 
             shape=20, col = "green") + 
  geom_point(aes(x=2, y=median(ratio_B[101:200,2])), shape=20, col = "green") +
  geom_hline(yintercept=1, lty=2) +
  ylim(0,2)+
  annotate("text", label=" Scenario B - Estimated CPDAG", x=2, y=1.7, size=4) +
  labs(x="", y="RMSE")




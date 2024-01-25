# This file contains the simulation study to compare the bayesian estimator based 
# on the posterior mean with the estimated average causal effects in Optimal IDA and 
# in the causal effect computed as f(L)- path (the methodology developed in this thesis)

source("auxiliary_functions.R")
source("final_function.R")

library(progress)
library(pcalg) 
library(MASS) 
library(parallel) 
library(abind) 
library(ggplot2) 
library(igraph)
library(BCDAG)


sim2 <- function(n, p, truecovm, x, y, x.lab, y.lab, d) {
  
  # generate data according to trueDAG
  dat <- mvtnorm::rmvnorm(n=n, mean=rep(0, p), sigma= unname(truecovm))
  colnames(dat) <- as.character(1:p)
  
  # estimate CPDAG
  score = new("GaussL0penObsScore", dat)
  estCPDAG = ges(score)
  estCPDAGamat = 1*wgtMatrix(estCPDAG$essgraph)
  dimnames(estCPDAGamat) = list(1:p, 1:p)
  
  # estimate the possible causal effects from the estimated CPDAG
  
  estOptimal  = ida(x.pos=x, y.pos=y, mcov=cov(dat), graphEst=t(estCPDAGamat), 
                    method="optimal")
  estPath     = causal_effect_path(p = p, CPDAGamat = t(estCPDAGamat),
                                   lab = colnames(dat), X = dat, x.pos = x, y.pos = y, 
                                   x.lab = x.lab, y.lab =y.lab)
  estBayesian = causal_effect_bayesian(X = dat, n = n, x.pos = x, y.pos = y, 
                                       p_edge = (d/p), S = 15000, b = 500)
  
  res = c(mean(estOptimal), mean(estPath), mean(estBayesian))
  return(res)
}


sim <- function(p, d, n){
  
  nonnull <- FALSE
  while (!nonnull) {
    
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
    
    
      # determine true smallest possible effect in CPDAG
    true <- ida(x, y, trueCov(trueDAG, back.compatible=TRUE), t(DAGamat), 
                  type="pdag")
      
    if (true>(10^(-7))) {nonnull <- TRUE}
    
  }
  
  res.def <- replicate(20, sim2(n , p, truecovm, x, y, x.lab, y.lab, d))
  
  MSE =  apply(res.def, 1, function (x) mean((x-true)^2) )
  names(MSE) <- c("MSE avg Opt", "MSE avg path", "MSE posterior mean")
  
  return(MSE)
  
}

# p = 10, d = 2, n= 500
simul.1 <- replicate(30, {
  sim(p=10, d=2, n=500)
})


write.csv(simul.1, "simul_1_bayesian.csv", row.names = F)


# Simul 1 - plot of MSE distribution of the thre methods under comparison
setwd("/Users/laura/Desktop/TESI/codes/results")
out = read.csv("simul_1_bayesian.csv")
median1 = apply(out, 1, function(x) median(x))
geomean1 = apply(out, 1, function(x) exp(mean(log(x))))
apply(out,1,mean)
ggdata_true <- data.frame( x=rep(c("MSE avg Opt", "MSE avg path", "MSE posterior mean"), each=30),
                           y=c(unlist(out[1, ]), unlist(out[2, ]),
                               unlist(out[3,])) )

ggdata_true$x <- factor(ggdata_true$x, 
                        levels=c("MSE avg Opt", "MSE avg path", "MSE posterior mean"))


ggplot(data=ggdata_true, aes(x=x, y=y)) +
  geom_boxplot() + 
  geom_point(aes(x=1, y=geomean1[1]), shape=20, col = "red") + 
  geom_point(aes(x=2, y=geomean1[2]), shape=20, col = "red") +
  geom_point(aes(x=3, y=geomean1[3]), shape=20, col = "red") +
  ylim(0, 0.1) +
  annotate("text", label="", x=2, y=4, size=6) +
  labs(x="", y="MSE of avg causal effect w.r.t. true causal effect") 

# Although we wouldn't say that there is a significant difference in the distribution of the MSE 
# associated with the three different methods, the bayesain estimator still exhibits 
# superior performances. This is probably due to its capability of handling model uncertainty

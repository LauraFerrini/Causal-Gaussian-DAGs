# Final function to find the causal effect from an estimated cpdag
# via our methodology 

library(purrr)
library(parallel)
library(BCDAG)

causal_effect_path <- function(p, CPDAGamat, lab, X, x.pos, y.pos, x.lab, y.lab) {
  
  out = MAXPDAGfromCPDAG_local(x.pos, y.pos, t(CPDAGamat))
  
  if (is.null(out)==FALSE){
    
    t.out = lapply(out, function(x) t(x)) # these are adjacency matrices of the dags encoded in our format

    undir_dags = mclapply(t.out, function(x) x * t(x), mc.cores = 5)
    dir_dags = Map("-", t.out, undir_dags)
    
    to_be_added = mclapply(undir_dags, function(x) {
      x[upper.tri(x)] = 0
      x}, mc.cores = 5)
  
    dags = Map("+", dir_dags, to_be_added)
    dags.lt = mclapply(dags, make_adj_lower_tri, p, lab, mc.cores = 5)
    CE = c()
    for (i in 1:length(dags.lt)){
      lab_i = colnames(dags.lt[[i]])
      
      CE[i] = find_CE_B(dags.lt[[i]], X[,lab_i], x.lab, y.lab)
      
    }
  }
  else{
    
    CE = 0
  }
  return(CE)
}



##########################
### BAYESIAN INFERENCE ###
##########################
to_list = function(x, along=length(dim(x))) {
  
  apply(x, along, identity, simplify=F)
}

causal_effect_bayesian <- function(X, p = ncol(X), n, S = 1000, b = 200, x.pos, y.pos, p_edge, 
                                   d = d) {
  
  out = learn_DAG(S = S, b = b, data = X, a = ncol(X), U = diag(1, ncol(X)),
                  w = p_edge)
  
  L_list = to_list(out$L)
  D_list = to_list(out$D)
  
  post_causal = unlist(Map(causaleffect, L_list, D_list, targets = x.pos, response = y.pos))
  
  return(post_causal)
}


# This code uses function implemented in the source code of 'pcalg' R package 
# to implement the causal effect estimation from the weight matrix L taking into account
# model uncertainty as well.



##########################
### DAGfromCPDAG_local ##
##########################

# Given the adjancency matrix of a cpdag (output of some structure learning alg)
# the algorithm proceeds semi-locally: by finding all the possible combinations
# of parents of the intervened node, then one can use meeks orientation rules to find the 
#  MAXpdags
# Comes from a readaptation of function optimal.est in github source code of 'pcalg'

MAXPDAGfromCPDAG_local <- function (x.pos, y.pos, amat.cpdag){
  
  
  amat.cpdag[which(amat.cpdag != 0)] <- 1  ##just in case make all the non-zero's 1
  amat.undir <- amat.cpdag * t(amat.cpdag) ## amat.undir - all undirected edges
  ## the undirected edge i - j has amat[i,j] = amat[j,i] =1
  
  amat.dir <- amat.cpdag - amat.undir  
  ## amat.dir - all directed edges
  pasets.dir <- which(amat.dir[x.pos,] != 0)   ## find all parents of x.pos in the PDAG
  
  ## sibs will be a vector containing all undirected edges connected with x.pos
  ## if for example: i is in x.pos and i-j is in G
  ## then add j;i to sibs
  sibs <- c()
  
  tmp.sib <- which(amat.undir[x.pos,]!=0)  ## tmp.sib contains all siblings of x.pos[i]
  
  if (length(tmp.sib)!=0){  ## if x.pos[i] has a sibling, then add all those sibling edges to sibs
    for (j in 1:length(tmp.sib)){
      sibs <- c(sibs,paste(tmp.sib[j],x.pos, sep = ";")) ## sibs is a vector of type a;b c;d
      ## where b and d are in x.pos and b-a d-c are in G
      ## note that if a,b are in x.pos and a-b is in G
      ## then both a;b and b;a are in sibs
    }
  }
  
  ## if sibs is empty, AND ch(X) = 0, then there is no possible path from x to y. return causal effect = to 0 
  
  if (length(sibs)==0 & length(which(amat.dir[,x.pos]!=0))==0){
    #cat("No possibly directed path from X to Y")
    #cat("the causal effect is 0")
    return(NULL)
  }
  
  if (length(sibs)==0){
    #cat("X has no siblings that need to be directed!")
    #cat("Use your input dag")
    return(list(amat.cpdag))
  } 
  
  size <- length(sibs)
  adj_matrices = list()
  m = 1
  
  ## first check the no additional parents option
  ## meaning that it is possible to orient all sibling edges out of x.pos
  toAdd <- mymakeBgKnowledge(sibs)
  
  
  amat.amenable <- myaddBgKnowledge(amat.cpdag, x = toAdd$y, y = toAdd$x,checkInput = FALSE)
  
  if (!is.null(amat.amenable)){
     
    adj_matrices[[m]] = amat.amenable
    m = m+ 1
    }
  
    
    
    ## now for all subsets of possible parents that is
    ## all possible subsets of sibs union pasets.dir:
    ## 1. check if they are allowed
    ## 2. orient all the possible edges by means of meeks orientation rules 
    
    
    for (k in 1:size){
      
      
      possPa <- combn(sibs,k) 
      
      ## form all combinations of sibs of size k
      for (r in 1:length(possPa[1,])){
        #r takes value from 1 to the number of possPA of size k
        s <- possPa[,r]   ## get one subset of sibs which we will try to orient into x.pos
        toAdd1 <- mymakeBgKnowledge(s)  ## transform it from sibs format: a;b c;d into bgKnowledge format that is
        ## into a data.frame with x=c(a,c) y=c(b,d) so that a -> b, c-> d is bgKnowledge
        
        sbar <- setdiff(sibs,s)       ## the complement of our subset s should be oriented out of x.pos
        toAdd2 <- mymakeBgKnowledge(sbar)
        
        ## the following 2 lines define the bg Knowledge that we try to add:
        #ie we orient possPA_x -> x, and the remaining siblings sbar as x-> sbar
        addFromThis <- c(toAdd1$x,toAdd2$y)
        addToThis <- c(toAdd1$y,toAdd2$x)
        
        
        ## only try to add this bg knowledge if its consistent within itself
        ## meaning if it does not contain contradictory edges
        ## for example addFromThis =c(1,2) addToThis = c(2,1)
        
        check2 <- FALSE  ## if check is true it will indicate that the background knowledge contradicts itself
        for (i in 1: length(addFromThis)){
          if (addToThis[i] %in% addFromThis[which(addToThis == addFromThis[i])]){# %in% addToThis)]){
            check2 <- TRUE
          }
        }
        if (!check2){  ##only try to add background knowledge that does not contradict itself
          amat.amenable <- myaddBgKnowledge(amat.cpdag, x=addFromThis,y=addToThis)
          if (!is.null(amat.amenable)){
            
            adj_matrices[[m]] = amat.amenable
            m = m +1
          }
        }
      }
    }
 
  return(adj_matrices)
}



#######################
### makeBgKnowledge ###
#######################

## the following functions takes character vector s of type s=c("1;2","3;4")
## and transforms it into a data frame containing vectors x=c(1,3) and y=c(2,4)
## so that one can add x -> y (that is 1 -> 2, 3 -> 4) as bg knowledge easier
## it is used in function MAXPDAGfromCPDAG to transform the siblings of x

mymakeBgKnowledge <- function(s)
{
  x <- y <- c()
  
  if (length(s)==0){   ## if s is empty, return empty vectors
    df.final <- data.frame(x=x, y=y)
    return(df.final)
  } else {            ## otherwise transform the charactor vector s into a
    ## bg knowledge data frame
    for (i in 1:length(s))
    {
      addFromTo <- as.numeric(unlist(strsplit(x = s[i],split = ";")))
      x <- c(x,addFromTo[1])
      y <- c(y,addFromTo[2])
    }
    df.final <- data.frame(x=x, y=y)
    return(df.final)
  }
}




#######################
### addBgKnowledge ###
######################


## This function is supposed to add the bg knowledge x -> y to the pdag in gInput
## x and y should be vectors of node positions. it checks for validity of the orientation
## by means of the Meeks orientation rules.


## gInput can be either the graphNEL object or adjacency matrix
## same encoding as in amat.cpdag above m[i,j]=0, m[j,i]=1 <=> i->j

myaddBgKnowledge <- function(gInput, x=c(), y=c(), verbose=FALSE, checkInput = TRUE){
  
  res <- gInput
  ##new line
  if (!is.matrix(gInput)) { #ie you have provided to the function a graphNEL obj
    
    if (numEdges(gInput)>0) {
      g <- t(as(gInput,"matrix")) ## g_ji if i->j # Converts the graphNEL into an adj matrix
      p <- as.numeric(dim(g)[1])
    }
    else {
      if (verbose) cat("Invalid or empty PDAG! This function only accepts graphNEL or adj mat\n")
      return(NULL)
    }
  }
  else {
    ##for me its easier to use an adjacency matrix
    if (length(res[1,])>0) {
      g <- res
      p <- length(g[1,])
    } else {
      if (verbose) cat("Invalid or empty PDAG! This function only accepts graphNEL or adj mat\n")
      return(NULL)
    }
  }
  pdag <- g
  #lab <- dimnames(g)[[1]]
  ## check if input is valid pdag
 
  ## real code starts from here
  ## previously was just dealing w gInput type
  ##how many new orientations -> k
  k <- length(x)
  i <- 1
  ## orient edge by edge
  ## after adding one orientation complete the orientation rules
  
  ##NEW if there are no edges to orient!!
  if (k ==0){
    
    pdag <- t(myapplyOrientationRules(t(pdag), verbose))
    if (!is.matrix(res)) {
      res <- as(t(pdag),"graphNEL")
    }
    else {
      res <- pdag
    }
    return(res)
  }
  
  
  while (i <=k){
    
    ## for each new edge
    ##find the nodes by the labels
    from <- x[i]
    to <- y[i]
    
    ##add the orientation if the edge in the current pdag is undirected
    ##and complete the orientation rules
    #if(length((pdag[from, to] == 1) & (pdag[to, from] == 1))!=0){
    if ((pdag[from,to] ==1) & (pdag[to,from] ==1))  {
      pdag[from,to] <- 0
      if (verbose) cat("Added orientation",x[i],"->",y[i],"to the PDAG.\n")
      pdag <- t(myapplyOrientationRules(t(pdag),verbose))  ##uses different matrix encoding
    #} 
    }else {
      ## the orientation we want to add conflicts with the current pdag
      ##either the opposite orientation is present in the current pdag
      ##or there is no edge between these two nodes in the pdag at all
      if ((pdag[from,to] ==1) & (pdag[to,from] ==0)){
        
        if (verbose) cat("Invalid bg knowledge! Cannot add orientation ",x[i],"->",y[i]," because",y[i],"->",x[i],"is already in the PDAG. \n")
        return(NULL)
      }
      
      if ((pdag[from,to] ==0) & (pdag[to,from] ==0)){
        if (verbose) cat("Invalid bg knowledge! Cannot add orientation",x[i],"->",y[i]," because there is no edge between",x[i],"and",y[i],"in the PDAG. \n")
        return(NULL)
      }
      
    }
    i <- i+1
  }
  
  if (!is.matrix(res))
  {
    res <- as(t(pdag),"graphNEL")
  } else {
    res <- pdag
  }
  return(res)
}



###############################
### Meeks orientation rules ###
###############################

myapplyOrientationRules <- function(gInput, verbose=FALSE) {
  res <- gInput
  ##new line
  if (!is.matrix(gInput))
  {
    if (numEdges(gInput)>0) {
      g <- as(gInput,"matrix") ## g_ij if i->j
      p <- as.numeric(dim(g)[1])
      pdag <- g
      ind <- which(g==1,arr.ind=TRUE)
    } else {
      cat("Invalid or empty PDAG! This function only accepts graphNEL or adj mat!\n")
      return(NULL)
    }
  } else {
    ##also new, for me its easier to use an adjacency matrix
    if (length(res[1,])>0){
      g <- res
      p <- length(g[1,])
      pdag <- g
      ind <- which(g==1,arr.ind=TRUE)
    } else {
      cat("Invalid or empty PDAG! This function only accepts graphNEL or adj mat!\n")
      return(NULL)
    }
  }
  ## Convert to complete pattern: use rules by Pearl/Meek also someone else Verma?
  old_pdag <- matrix(0, p,p)
  while (!all(old_pdag == pdag)) {
    old_pdag <- pdag
    ## rule 1
    ind <- which((pdag==1 & t(pdag)==0), arr.ind=TRUE) ## a -> b
    for (i in seq_len(nrow(ind))) {
      a <- ind[i,1]
      b <- ind[i,2]
      indC <- which( (pdag[b,]==1 & pdag[,b]==1) & (pdag[a,]==0 & pdag[,a]==0))
      if (length(indC)>0) {
        pdag[b,indC] <- 1
        pdag[indC,b] <- 0
        if (verbose)
          cat("\nRule 1:",a,"->",b," and ",b,"-",indC,
              " where ",a," and ",indC," not connected: ",b,"->",indC,"\n")
      }
    }
    ## x11()
    ## plot(as(pdag,"graphNEL"), main="After Rule1")
    
    ## rule 2
    ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a -> b
    for (i in seq_len(nrow(ind))) {
      a <- ind[i,1]
      b <- ind[i,2]
      indC <- which( (pdag[a,]==1 & pdag[,a]==0) & (pdag[,b]==1 & pdag[b,]==0))
      if (length(indC)>0) {
        pdag[a,b] <- 1
        pdag[b,a] <- 0
        if (verbose) cat("\nRule 2: Kette ",a,"->",indC,"->",
                         b,":",a,"->",b,"\n")
      }
    }
    ## x11()
    ## plot(as(pdag,"graphNEL"), main="After Rule2")
    
    ## rule 3
    ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
    for (i in seq_len(nrow(ind))) {
      a <- ind[i,1]
      b <- ind[i,2]
      indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==1 & pdag[b,]==0))
      if (length(indC)>=2) {
        ## cat("R3: indC = ",indC,"\n")
        g2 <- pdag[indC,indC]
        ## print(g2)
        if (length(g2)<=1) {
          g2 <- 0
        } else {
          diag(g2) <- rep(1,length(indC)) ## no self reference
        }
        ##Bugfix: changed  #####
        ## if (any(g2 == 0)) {
        g3 <- g2 + t(g2)
        if (any(g3 == 0)) {
          ###### end of change
          ## if (any(g2==0)) { ## if two nodes in g2 are not connected
          pdag[a,b] <- 1
          pdag[b,a] <- 0
          if (verbose) cat("\nRule 3:",a,"->",b,"\n")
        }
      }
    }
    ## x11()
    ## plot(as(pdag,"graphNEL"), main="After Rule3")
    
    ## rule 4
    ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
    if (length(ind)>0) {
      for (i in seq_len(nrow(ind))) {
        a <- ind[i,1]
        b <- ind[i,2]
        indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==0 & pdag[b,]==0))
        l.indC <- length(indC)
        if (l.indC>0) {
          found <- FALSE
          ic <- 0
          while(!found & (ic < l.indC)) {
            ic <- ic + 1
            c <- indC[ic]
            indD <- which( (pdag[c,]==1 & pdag[,c]==0) & (pdag[,b]==1 & pdag[b,]==0))
            if (length(indD)>0) {
              found <- TRUE
              pdag[b,a] = 0
              if (verbose) cat("Rule 4 applied \n")
            }
          }
        }
      }
      
      
    }
  }
  if (!is.matrix(res))
  {
    res <- as(pdag,"graphNEL")
  } else {
    res <- pdag
  }
  return(res)
}




#########################
## make_adj_lower_tri ##
########################

# The following function takes as input an adjacency matrix of a dag (encoded in our format)
# and transforms it in its lower triangular version

## nb: adj_matrix must be provided without labels

make_adj_lower_tri <- function(adj_matrix, p, lab) {
  # Initialization
  adj_lowert = matrix(0, p, p)
  my_list = list()
  # create a list such that each element my_list[[i]] contains the parents of node i
  for (i in 1:p){
    my_list[[i]] = which(adj_matrix[,i]!=0)
  }
  
  npa = colSums(adj_matrix)
  # Initialize the vector containing the "new" column order so that in the first position there is the column having
  # 0 parents
  
  order_def = unname(c(which(npa==0)))
 
  for (m in order_def){
    my_list[[m]] = 0
  }
  
  k = length(order_def)
  
  while (k<p) {
    for (i in 1:p){
      
      if (length(setdiff(my_list[[i]], order_def[1:k]))==0){
        
        order_def[k+1] = i
        k = k+1
        my_list[[i]] = 0 
        break
      }
    }
    
  }
  
  adj_lowert[1:p, 1:p] = (adj_matrix[rev(order_def), rev(order_def)]) # columns and row exchange
  dimnames(adj_lowert) = list( lab[rev(order_def)], lab[rev(order_def)]) # assign new labels
  return(adj_lowert)
}



#################
## find_CE_B ####
#################

# INPUT: X: data matrix, adj_DAG: matrix in lower triangular format (with dimnames),
# x.lab, y:lab: x, and y labels

find_CE_B = function(adj_DAG, X, x.lab, y.lab){
  lab   = rownames(adj_DAG)
  x.pos = which(lab==x.lab)
  y.pos = which(lab==y.lab)
    
  Mu_mle   = colMeans(X)
  Sigma_mle = matrix(0, dim(X)[2], dim(X)[2])
  B_hat = matrix(0, dim(X)[2], dim(X)[2])
    for (i in 1:(dim(X)[1])){
      Sigma_mle = Sigma_mle + (X[i,]-Mu_mle)%*%t(X[i,]-Mu_mle)
      
    }
  
  Sigma_mle = Sigma_mle/(dim(X)[1])
  
  for (i in 1:(dim(X)[2])){
    pa_i = which(adj_DAG[,i] == 1)
    if(length(pa_i!=0)){
      B_hat[pa_i, i] = (Sigma_mle[i, pa_i] %*% solve(Sigma_mle[pa_i,pa_i])) # notice that B_hat = I-L
    }
  
  }
  ce_est = B_hat[x.pos,y.pos]
  j = 1
  
  B_hatj = B_hat
  
  while(j <= x.pos-1){
  
      B_hatj = B_hatj  %*% B_hat
      ce_est = ce_est + B_hatj[x.pos, y.pos]
      
      j = j+1
  }
  return(ce_est)
  
}




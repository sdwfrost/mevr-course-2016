slidingWindowAlignment <- function(s,l,k){
  seqlen <- dim(s)[[2]]
  startpos <- 1
  endpos <- 1+l-1
  result <- list()
  i <- 1
  while(endpos<=seqlen){
    result[[i]] <- s[,startpos:endpos]
    i <- i+1
    startpos <- startpos+k
    endpos <- endpos+k
  }
  result
}

sbtest <- function(s,minl,model){
  seqlen <- dim(s)[[2]]
  breakpos <- minl+1
  bp <- c() 
  symdiff <- c()
  bsdiff <- c()
  pathdiff <- c()
  quadpathdiff <- c()
  rfdist <- c()
  while(breakpos<seqlen-minl){
    # print(breakpos)
    lhs <- s[,1:breakpos]
    rhs <- s[,(breakpos+1):seqlen]
    lhs.dist <- dist.dna(lhs,model=model,as.matrix=TRUE)
    rhs.dist <- dist.dna(rhs,model=model,as.matrix=TRUE)
    # need to check why nj doesn't work here
    lhs.nj <- njs(lhs.dist)
    rhs.nj <- njs(rhs.dist)
    d1 <- treedist(lhs.nj,rhs.nj)
    bp <- c(bp,breakpos)
    symdiff <- c(symdiff,d1[1])
    bsdiff <- c(bsdiff,d1[2])
    pathdiff <- c(pathdiff,d1[3])
    quadpathdiff <- c(quadpathdiff,d1[4])
    breakpos <- breakpos+1
  }
  data.frame(breakpoint=bp,symdiff,bsdiff,pathdiff,quadpathdiff)
}

disttree <- function(l){
  numtrees <- length(l)
  result <- list()
  result[["symmetric.difference"]] <- matrix(0,numtrees,numtrees)
  result[["branch.score.difference"]] <- matrix(0,numtrees,numtrees)
  result[["path.difference"]] <- matrix(0,numtrees,numtrees)
  result[["quadratic.path.difference"]] <- matrix(0,numtrees,numtrees)
  for(i in 1:(numtrees-1)){
    for(j in (i+1):numtrees){
      d <- treedist(l[[i]],l[[j]])
      for(k in 1:4){
        result[[k]][i,j] <- d[k]
        result[[k]][j,i] <- d[k]
      }
    }
  }
  result
}
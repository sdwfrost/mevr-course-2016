require(ape)
require(INLA)

branching.sampling.times <- function(phy){
  phy = new2old.phylo(phy)
  if (class(phy) != "phylo")
    stop("object \"phy\" is not of class \"phylo\"")
  tmp <- as.numeric(phy$edge)
  nb.tip <- max(tmp)
  nb.node <- -min(tmp)
  xx <- as.numeric(rep(NA, nb.tip + nb.node))
  names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
  xx["-1"] <- 0
  for (i in 2:length(xx)) {
    nod <- names(xx[i])
    ind <- which(phy$edge[, 2] == nod)
    base <- phy$edge[ind, 1]
    xx[i] <- xx[base] + phy$edge.length[ind]
  }
  depth <- max(xx)
  branching.sampling.times <- depth - xx
  return(branching.sampling.times)
}

heterochronous.gp.stat <- function(phy){
  b.s.times = branching.sampling.times(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  coal.interval = rep(NA, num.coal.events)
  coal.lineages = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  #unique.sampling.times = sort(unique(b.s.times[tip.ind]))
  sampling.times = sort((b.s.times[tip.ind]))
  for (i in 2:length(sampling.times)){
   if ((sampling.times[i]-sampling.times[i-1])<0.1){
     sampling.times[i]<-sampling.times[i-1]}
  }
  unique.sampling.times<-unique(sampling.times)
  sampled.lineages = NULL
  for (sample.time in unique.sampling.times){
   sampled.lineages = c(sampled.lineages,
    sum(sampling.times == sample.time))  
  }
return(list(coal.times=sorted.coal.times, sample.times = unique.sampling.times, sampled.lineages=sampled.lineages))  
}

calculate.skyride <- function(tree){
  ##Data prep to use INLA
  dd<-heterochronous.gp.stat(tree)
  s.time<-dd$sample.times
  n.sample<-dd$sampled.lineages
  n<-length(dd$coal.times)+1
  data<-matrix(0,nrow=n-1,ncol=2)
  data[,1]<-dd$coal.times
  s.time<-c(s.time,max(data[,1])+1)
  data[1,2]<-sum(n.sample[s.time<=data[1,1]])
  tt<-length(s.time[s.time<=data[1,1]])+1
  for (j in 2:nrow(data)){
    if (data[j,1]<s.time[tt]){
      data[j,2]<-data[j-1,2]-1
		}else{
			data[j,2]<-data[j-1,2]-1+sum(n.sample[s.time>data[j-1,1] & s.time<=data[j,1]])
		  tt<-length(s.time[s.time<=data[j,1]])+1
		}	
  }
  ###Bayesian Skyride -CGGP
  s<-unique(sort(c(data[,1],s.time[1:length(s.time)-1])))
  event1<-sort(c(data[,1],s.time[1:length(s.time)-1]),index.return=TRUE)$ix
  n<-nrow(data)+1
  l<-length(s)
  event<-rep(0,l)
  event[event1<n]<-1
  y<-diff(s)
  coal.factor<-rep(0,l-1)
  indicator<-rep(0,l-1)
  t<-rep(0,l-1)
  indicator<-cumsum(n.sample[s.time<data[1,1]])
  indicator<-c(indicator,indicator[length(indicator)]-1)
  ini<-length(indicator)+1
  for (k in ini:(l-1)){
    j<-data[data[,1]<s[k+1] & data[,1]>=s[k],2]
    if (length(j)==0){indicator[k]=indicator[k-1]+sum(n.sample[s.time<s[k+1] & s.time>=s[k]])}
    if (length(j)>0){indicator[k]<-j-1+sum(n.sample[s.time<s[k+1] & s.time>=s[k]])}
  }
  coal.factor<-indicator*(indicator-1)/2
  coal.times<-data[,1]
  u<-c(coal.times[1],diff(coal.times))
  label<-c(1,1+cumsum(event[1:(l-1)]))
  factor.matrix<-.5*diag(length(u))
  for (j in 2:length(u)){
    for (i in 1:(j-1)){
      factor.matrix[j,i]<-1
    }
  }
  time_aware_short<-factor.matrix%*%u
  time_aware<-time_aware_short[label]
  E <- y*coal.factor
  formula=y~-1+f(time,model="rw1",hyper = list(prec = list(param = c(.001, .001))),constr=FALSE)
  data_mod<-data.frame(y=event[-1],event=event[-1],time=time_aware[-1],E=log(E))
  mod5<-inla(formula,family="poisson",data=data_mod,offset=E,control.predictor=list(compute=TRUE))
  grid.points.scaled1<-dd$coal.times
  sr.median <- exp(-mod5$summary.random$time$"0.5quant")
  sr.lc <- exp(-mod5$summary.random$time$"0.025quant")
  sr.uc <- exp(-mod5$summary.random$time$"0.975quant")
  data.frame(time=grid.points.scaled1,sr.median=sr.median,sr.lc=sr.lc,sr.uc=sr.uc)
}


require(seqinr)

countstops <-
function(ch,sens="F",numcode=1){
  # Calculate translations
  if(sens=="F"|sens=="R"){
    f0 <- seqinr::translate(ch,sens=sens,numcode=numcode)
    f2 <- seqinr::translate(c("-",ch),sens=sens,numcode=numcode)
    f1 <- seqinr::translate(c("-","-",ch),sens=sens,numcode=numcode)
  }
  else if(sens=="B"){
    f0 <- seqinr::translate(ch,sens="F",numcode=numcode)
    f2 <- seqinr::translate(c("-",ch),sens="F",numcode=numcode)
    f1 <- seqinr::translate(c("-","-",ch),sens="F",numcode=numcode)
    r0 <- seqinr::translate(ch,sens="R",numcode=numcode)
    r2 <- seqinr::translate(c("-",ch),sens="R",numcode=numcode)
    r1 <- seqinr::translate(c("-","-",ch),sens="R",numcode=numcode)
  }
  f0stops <- sum(f0=="*")
  f1stops <- sum(f1=="*")
  f2stops <- sum(f2=="*")
  if(sens=="B"){
    r0stops <- sum(r0=="*")
    r1stops <- sum(r1=="*")
    r2stops <- sum(r2=="*")
    return(c(f0stops,f1stops,f2stops,r0stops,r1stops,r2stops))
  }
  else{
    return(c(f0stops,f1stops,f2stops))
  }
}


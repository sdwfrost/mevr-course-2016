require(XML)
require(magrittr)
require(annotate)

accToIds <- function(accs){
  myids <- entrez_fetch(db="nuccore",id=accs,rettype="uilist") %>%
    xmlParse  %>%
    xmlToList %>%
    unlist %>%
    unname
}

null.to.other <-
function(x,y=NA){
  if(is.null(x)){
    return(y)
  }
  else{
    return(x)
  }
}

read.gbc <- function(fn,verbose=TRUE){
  myxml <- xmlTreeParse(fn)
  xmltop <- xmlRoot(myxml)
  numseq <- length(xmltop)
  accession <- rep("",numseq)
  country <- rep(NA,numseq)
  collectiondate <- rep(NA,numseq)
  host <- rep(NA,numseq)
  for(i in 1:numseq){
    if(verbose){
      print(paste("Processing",i,"of",numseq))
    }
    accession[i] <-  unlist(lapply(getNodeSet(xmltop[[i]],"//INSDSeq/INSDSeq_primary-accession"),xmlValue))[1]
    country[i] <- null.to.other(unlist(lapply(getNodeSet(xmltop[[i]],"//INSDSeq/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier[INSDQualifier_name='country']/INSDQualifier_value"),xmlValue))[1],NA)
    collectiondate[i] <- null.to.other(unlist(lapply(getNodeSet(xmltop[[i]],"//INSDSeq/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier[INSDQualifier_name='collection_date']/INSDQualifier_value"),xmlValue))[1],NA)
    host[i] <- null.to.other(unlist(lapply(getNodeSet(xmltop[[i]],"//INSDSeq/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier[INSDQualifier_name='host']/INSDQualifier_value"),xmlValue))[1],NA)
  }
    return(data.frame(accession,country,collectiondate,host))
}

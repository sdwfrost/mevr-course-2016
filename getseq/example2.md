# Getting annotations

Sequences in Genbank also have metadata associated with them. I have included a function, ```read.gbc``` in ```utils.R``` that extracts the accession, country, collection date and host from sequence data.


```r
source("utils.R")
```

This function is rather slow at the moment, but works on datasets of less than 500 sequences.


```r
annotations <- read.gbc("denv.xml")
```

```
## [1] "Processing 1 of 116"
## [1] "Processing 2 of 116"
## [1] "Processing 3 of 116"
## [1] "Processing 4 of 116"
## [1] "Processing 5 of 116"
## [1] "Processing 6 of 116"
## [1] "Processing 7 of 116"
## [1] "Processing 8 of 116"
## [1] "Processing 9 of 116"
## [1] "Processing 10 of 116"
## [1] "Processing 11 of 116"
## [1] "Processing 12 of 116"
## [1] "Processing 13 of 116"
## [1] "Processing 14 of 116"
## [1] "Processing 15 of 116"
## [1] "Processing 16 of 116"
## [1] "Processing 17 of 116"
## [1] "Processing 18 of 116"
## [1] "Processing 19 of 116"
## [1] "Processing 20 of 116"
## [1] "Processing 21 of 116"
## [1] "Processing 22 of 116"
## [1] "Processing 23 of 116"
## [1] "Processing 24 of 116"
## [1] "Processing 25 of 116"
## [1] "Processing 26 of 116"
## [1] "Processing 27 of 116"
## [1] "Processing 28 of 116"
## [1] "Processing 29 of 116"
## [1] "Processing 30 of 116"
## [1] "Processing 31 of 116"
## [1] "Processing 32 of 116"
## [1] "Processing 33 of 116"
## [1] "Processing 34 of 116"
## [1] "Processing 35 of 116"
## [1] "Processing 36 of 116"
## [1] "Processing 37 of 116"
## [1] "Processing 38 of 116"
## [1] "Processing 39 of 116"
## [1] "Processing 40 of 116"
## [1] "Processing 41 of 116"
## [1] "Processing 42 of 116"
## [1] "Processing 43 of 116"
## [1] "Processing 44 of 116"
## [1] "Processing 45 of 116"
## [1] "Processing 46 of 116"
## [1] "Processing 47 of 116"
## [1] "Processing 48 of 116"
## [1] "Processing 49 of 116"
## [1] "Processing 50 of 116"
## [1] "Processing 51 of 116"
## [1] "Processing 52 of 116"
## [1] "Processing 53 of 116"
## [1] "Processing 54 of 116"
## [1] "Processing 55 of 116"
## [1] "Processing 56 of 116"
## [1] "Processing 57 of 116"
## [1] "Processing 58 of 116"
## [1] "Processing 59 of 116"
## [1] "Processing 60 of 116"
## [1] "Processing 61 of 116"
## [1] "Processing 62 of 116"
## [1] "Processing 63 of 116"
## [1] "Processing 64 of 116"
## [1] "Processing 65 of 116"
## [1] "Processing 66 of 116"
## [1] "Processing 67 of 116"
## [1] "Processing 68 of 116"
## [1] "Processing 69 of 116"
## [1] "Processing 70 of 116"
## [1] "Processing 71 of 116"
## [1] "Processing 72 of 116"
## [1] "Processing 73 of 116"
## [1] "Processing 74 of 116"
## [1] "Processing 75 of 116"
## [1] "Processing 76 of 116"
## [1] "Processing 77 of 116"
## [1] "Processing 78 of 116"
## [1] "Processing 79 of 116"
## [1] "Processing 80 of 116"
## [1] "Processing 81 of 116"
## [1] "Processing 82 of 116"
## [1] "Processing 83 of 116"
## [1] "Processing 84 of 116"
## [1] "Processing 85 of 116"
## [1] "Processing 86 of 116"
## [1] "Processing 87 of 116"
## [1] "Processing 88 of 116"
## [1] "Processing 89 of 116"
## [1] "Processing 90 of 116"
## [1] "Processing 91 of 116"
## [1] "Processing 92 of 116"
## [1] "Processing 93 of 116"
## [1] "Processing 94 of 116"
## [1] "Processing 95 of 116"
## [1] "Processing 96 of 116"
## [1] "Processing 97 of 116"
## [1] "Processing 98 of 116"
## [1] "Processing 99 of 116"
## [1] "Processing 100 of 116"
## [1] "Processing 101 of 116"
## [1] "Processing 102 of 116"
## [1] "Processing 103 of 116"
## [1] "Processing 104 of 116"
## [1] "Processing 105 of 116"
## [1] "Processing 106 of 116"
## [1] "Processing 107 of 116"
## [1] "Processing 108 of 116"
## [1] "Processing 109 of 116"
## [1] "Processing 110 of 116"
## [1] "Processing 111 of 116"
## [1] "Processing 112 of 116"
## [1] "Processing 113 of 116"
## [1] "Processing 114 of 116"
## [1] "Processing 115 of 116"
## [1] "Processing 116 of 116"
```

This function returns a dataframe of the annotations.


```r
head(annotations)
```

```
##   accession              country collectiondate         host
## 1  KC333651                China           2012 Homo sapiens
## 2  JX024758            Singapore       Dec-2010 Homo sapiens
## 3  JX024757            Singapore       Dec-2010 Homo sapiens
## 4  JQ513345     Brazil: Bahia-BA    18-Mar-2011 Homo sapiens
## 5  JQ513343    Brazil: Manaus-AM    29-Jan-2011 Homo sapiens
## 6  JQ513341 Brazil: Boa Vista-RR    21-Nov-2010 Homo sapiens
```


For those of you interested, here is the code.


```r
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
```

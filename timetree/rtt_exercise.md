# Root-to-tip exercise






```r
library(ape)
library(stringr)
library(phangorn)
source("chronos.R")
```


```r
seq.fn <- "H5N1.fas"
```


```r
stub <- strsplit(seq.fn,".",fixed=TRUE)[1]
myseqs <- read.dna(seq.fn,format="fasta",as.matrix=FALSE)
```

ML tree.


```r
njtree  <- nj(dist.dna(myseqs,model="TN93"))
myseqs.phydat <- as.phyDat(myseqs)
myseqs.gtrig <- pml(njtree,myseqs.phydat,model="GTR+I+G",k=4)
```

```
## Warning in pml(njtree, myseqs.phydat, model = "GTR+I+G", k = 4): negative
## edges length changed to 0!
```

```r
myseqs.gtrig <- optim.pml(myseqs.gtrig,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE,optRate=FALSE)
```

```
## optimize edge weights:  -6050.051 --> -6038.844 
## optimize base frequencies:  -6038.844 --> -5986.777 
## optimize rate matrix:  -5986.777 --> -5728.631 
## optimize invariant sites:  -5728.631 --> -5717.886 
## optimize shape parameter:  -5717.886 --> -5717.886 
## optimize edge weights:  -5717.886 --> -5717.684 
## optimize topology:  -5717.684 --> -5704.969 
## optimize topology:  -5704.969 --> -5704.969 
## 4 
## optimize base frequencies:  -5704.969 --> -5704.725 
## optimize rate matrix:  -5704.725 --> -5704.684 
## optimize invariant sites:  -5704.684 --> -5704.657 
## optimize shape parameter:  -5704.657 --> -5704.657 
## optimize edge weights:  -5704.657 --> -5704.657 
## optimize topology:  -5704.657 --> -5704.657 
## 0 
## optimize base frequencies:  -5704.657 --> -5704.655 
## optimize rate matrix:  -5704.655 --> -5704.655 
## optimize invariant sites:  -5704.655 --> -5704.655 
## optimize shape parameter:  -5704.655 --> -5704.655 
## optimize edge weights:  -5704.655 --> -5704.655 
## optimize base frequencies:  -5704.655 --> -5704.655 
## optimize rate matrix:  -5704.655 --> -5704.655 
## optimize invariant sites:  -5704.655 --> -5704.655 
## optimize shape parameter:  -5704.655 --> -5704.654 
## optimize edge weights:  -5704.654 --> -5704.654 
## optimize base frequencies:  -5704.654 --> -5704.654 
## optimize rate matrix:  -5704.654 --> -5704.654 
## optimize invariant sites:  -5704.654 --> -5704.654 
## optimize shape parameter:  -5704.654 --> -5704.654 
## optimize edge weights:  -5704.654 --> -5704.654 
## optimize base frequencies:  -5704.654 --> -5704.654 
## optimize rate matrix:  -5704.654 --> -5704.654 
## optimize invariant sites:  -5704.654 --> -5704.654 
## optimize shape parameter:  -5704.654 --> -5704.653 
## optimize edge weights:  -5704.653 --> -5704.653 
## optimize base frequencies:  -5704.653 --> -5704.653 
## optimize rate matrix:  -5704.653 --> -5704.653 
## optimize invariant sites:  -5704.653 --> -5704.653 
## optimize shape parameter:  -5704.653 --> -5704.653 
## optimize edge weights:  -5704.653 --> -5704.653 
## optimize base frequencies:  -5704.653 --> -5704.653 
## optimize rate matrix:  -5704.653 --> -5704.653 
## optimize invariant sites:  -5704.653 --> -5704.653 
## optimize shape parameter:  -5704.653 --> -5704.653 
## optimize edge weights:  -5704.653 --> -5704.653 
## optimize base frequencies:  -5704.653 --> -5704.653 
## optimize rate matrix:  -5704.653 --> -5704.653 
## optimize invariant sites:  -5704.653 --> -5704.653 
## optimize shape parameter:  -5704.653 --> -5704.652 
## optimize edge weights:  -5704.652 --> -5704.652 
## optimize base frequencies:  -5704.652 --> -5704.652 
## optimize rate matrix:  -5704.652 --> -5704.652 
## optimize invariant sites:  -5704.652 --> -5704.652 
## optimize shape parameter:  -5704.652 --> -5704.652 
## optimize edge weights:  -5704.652 --> -5704.652
```

```r
tr <- myseqs.gtrig$tree
```


```r
tipnames <- tr$tip.label
tipdates <- as.integer(str_sub(tipnames,-4))
```


```r
tr.rtt <- rtt(tr,tipdates)
tr.rtt$edge.length[tr.rtt$edge.length<0.00000001] <- 0.00000001
tr.rtt.tipnames <- tr.rtt$tip.label
tr.rtt.tipdates <- as.integer(str_sub(tr.rtt.tipnames,-4))
```



```r
rootdistance <- distRoot(tr.rtt)
pathlm <- lm(rootdistance~tr.rtt.tipdates)
rate <- coef(pathlm)[2]
rate
```

```
## tr.rtt.tipdates 
##    -0.001775018
```

```r
tmrca <- unname(-coef(pathlm)[1]/coef(pathlm)[2])
tmrca
```

```
## [1] 2024.524
```


```r
max.time <- max(tipdates)
ncat <- 1
strict.clock.ctrl <- chronos.control(nb.rate.cat=as.integer(ncat))
calibrating.values <- makeChronosCalib(tr)
calibrating.values$age.min <- max.time - root.time
calibrating.values$age.max <- max.time - root.time
# pins the tips to sampling years
calibrating.values <- rbind(calibrating.values,
                            data.frame(node=seq(1,length(td)),
                                       age.min=max.time - td,
                                       age.max=max.time - td,
                                       soft.bounds=FALSE))
dated.tree <- RLchronos(tr, 
                     lambda=1, 
                     model=sub.rate.model, 
                     calibration=calibrating.values,
                     control=strict.clock.ctrl)
```

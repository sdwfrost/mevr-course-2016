# Practice



## Load libraries


```r
library(ape)
library(stepwise)
library(magrittr)
library(kdetrees)
library(phangorn)
library(bios2mds)
```

```
## Loading required package: amap
## Loading required package: e1071
## Loading required package: scales
## Loading required package: cluster
## Loading required package: rgl
```

```r
source("recombination.R")
```

In addition to the libraries, I have also written some functions:

- ```slidingWindowAlignment```: generates sliding window alignments
- ```sbtest```: generates trees either side of a breakpoint, and calculates the distance between the trees
- ```disttree```: calculates distances between all trees in a list

## Set filename and outgroup


```r
seqfilename <- "CRF7.fas"
outgroup <- "J_SE7887"
```

## Read FASTA file


```r
seqdata <- read.dna(seqfilename,format="fasta",as.matrix=TRUE)
seqdata
```

```
## 10 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 7894 
## 
## Labels: C_C2220 CRF7_C54A A_92UG037 B_JRFL D_NDK F_93BR020 ...
## 
## Base composition:
##     a     c     g     t 
## 0.366 0.173 0.236 0.225
```

## Output data with 'clean' sequence names


```r
seqdata.stepwise <- seqdata
row.names(seqdata.stepwise) <- makeLabel(row.names(seqdata),len=9)
write.dna(seqdata.stepwise,file=paste(seqfilename,".stepwise",sep=""),format="interleaved",colsep="",nbcol=-1)
```

## Test for recombination using Maxchi


```r
seqdata.maxchi <- maxchi(paste(seqfilename,".stepwise",sep=""),breaks=seg.sites(seqdata),winHalfWidth=50,permReps=100)
summary(seqdata.maxchi)
```

```
## --------------------------------------------------
## There were 12 site-specific MaxChi statistics significant at the
## 10 percent level (90th percentile = 21.007, 95th percentile = 22.243):
## 
## Number Location  MaxChi   pairs
##      1      293  23.645*   (B_JRFL    :G_92NG083 )
##      2      294  21.583*   (B_JRFL    :G_92NG083 )
##      3      297  23.377*   (B_JRFL    :G_92NG083 )
##      4      300  25.253*   (B_JRFL    :G_92NG083 )
##      5      301  23.377*   (B_JRFL    :G_92NG083 )
##      6      387  21.236*   (B_JRFL    :G_92NG083 )
##      7     6022  23.188*   (CRF7_C54A :G_92NG083 )
##      8     6084  21.374*   (C_C2220   :G_92NG083 )
##      9     6085  23.377*   (C_C2220   :G_92NG083 )
##     10     6090  21.374*   (CRF7_C54A :G_92NG083 )
##     11     6091  23.188*   (CRF7_C54A :G_92NG083 )
##     12     6092  21.374*   (C_C2220   :G_92NG083 )
## --------------------------------------------------
## Notes - "Location" is the polymorphic site just before the proposed breakpoint.
##       - MaxChi statistics significant at the 5 percent level indicated by a *
```

## Test for recombination using phylogenetic profiling


```r
seqdata.phylpro <- phylpro(paste(seqfilename,".stepwise",sep=""),breaks=seg.sites(seqdata),winHalfWidth=100,permReps=100)
summary(seqdata.phylpro)
```

```
## Length  Class   Mode 
##      0   NULL   NULL
```

## Sliding window phylogeny

### Making a single tree

Making a tree with a sequence alignment using a distance based approach involves making a distance matrix, performing tree reconstruction, perhaps rooting the tree, and plotting it out.

The R library ```magrittr``` makes this more straightforward, by introducing an operator, ```%>%```, which sends the output of one command to the next. If the second command has more than one argument, we use the placeholder ```.``` instead. In this way, we can develop a mini-pipeline.


```r
dist.dna(seqdata,"TN93") %>% # Make distance matrix
  nj %>% # Neighbour joining
  root(.,outgroup=outgroup) %>% # Root with outgroup
  plot # Plot tree
```

```
## Error in if (newroot == ROOT) {: argument is of length zero
```

### Making a set of trees

As the sequences are stored in a matrix, it is straightforward to make a list of alignments, each of which is a window on the original alignment.

The following command generates sequence alignments 300 base pairs long, moving in steps of 10.


```r
seqdata.slide <- slidingWindowAlignment(seqdata,300,10)
length(seqdata.slide)
```

```
## [1] 760
```

This generates a list of alignments. In R, there is a command, ```lapply```, that applies a command to each element in a list. The following generates a list of distance matrices, then a list of neighbour joining trees.



```r
seqdata.slide.nj <- lapply(seqdata.slide,dist.dna,model="TN93",as.matrix=TRUE) %>%
  lapply(.,njs)
```

### Outlier detection

We now have a list of trees. How do we determine whether a single tree explains all the sub-alignments? One approach is to work out whether there are 'outlying trees'. The command ```kdetrees``` computes a distribution of trees, then determines whether there are trees in the 'tail' of the distribution.


```r
seqdata.slide.nj.kde <- kdetrees(seqdata.slide.nj,outgroup=outgroup)
```

Plotting out the output from this function will show the outlying trees, as a function of the index of the tree - each of which represents the tree from a slice of the original alignment.


```r
plot(seqdata.slide.nj.kde)
```

```
## Error in plot(seqdata.slide.nj.kde): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'seqdata.slide.nj.kde' not found
```

We can also plot out a histogram.


```r
hist(seqdata.slide.nj.kde)
```

```
## Error in hist(seqdata.slide.nj.kde): object 'seqdata.slide.nj.kde' not found
```

### Distances between trees

Another way is to work out the distances between trees.


```r
seqdata.slide.treedist <- disttree(seqdata.slide.nj)
```


```r
seqdata.slide.treemds <- mmds(seqdata.slide.treedist[[1]],pc=2)
mmds.2D.plot(seqdata.slide.treemds)
```

### Single breakpoint test

This function splits the alignment into two (at least 300 base pairs long), and calculates the distance between the trees for either side of the breakpoint.


```r
seqdata.sbt <- sbtest(seqdata,300,"TN93")
```

Now we can plot the distance between the trees as a function of the breakpoint.


```r
plot(symdiff~breakpoint,data=seqdata.sbt,type="s")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png) 

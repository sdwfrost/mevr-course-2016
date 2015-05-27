# Model choice exercise



Load the libraries.


```r
library(ape)
library(phangorn)
```

Load in the sequence data and the neighbour joining tree.


```r
myalignment <- read.dna("myseqs_aligned_pruned.fas",format="fasta",as.matrix=TRUE)
mytree  <- read.tree("myseqs_aligned_pruned_nj_tn93.tre")
```

Now compare different models.


```r
myseq.phydat <- as.phyDat(myseq) # Convert to a format modeltest understand
```

```
## Error in as.phyDat(myseq): object 'myseq' not found
```

```r
myseq.modeltest <- modelTest(myseq.phydat,tree=mytree,model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, k = 4, control = pml.control(epsilon = 1e-08, maxit = 3, trace = 1), multicore = FALSE)
```

```
## Error in modelTest(myseq.phydat, tree = mytree, model = c("JC", "F81", : object 'myseq.phydat' not found
```

```r
myseq.modeltest
```

```
## Error in eval(expr, envir, enclos): object 'myseq.modeltest' not found
```

Now obtain a maximum likelihood tree, starting from the neighbour joining tree. Although this can be done all in one command, using ```optim.pml```, this shows how you can get finer control, and in principle get a tree faster.


```r
myseq.gtrig <- pml(mytree,myseq.phydat,model="GTR+I+G",k=4)
```

```
## Warning in pml(mytree, myseq.phydat, model = "GTR+I+G", k = 4): negative
## edges length changed to 0!
```

```
## Error in pml(mytree, myseq.phydat, model = "GTR+I+G", k = 4): object 'myseq.phydat' not found
```

```r
# optimise
myseq.gtrig <- optim.pml(myseq.gtrig,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE,optRate=TRUE)
```

```
## Error in optim.pml(myseq.gtrig, optNni = TRUE, optBf = TRUE, optQ = TRUE, : object 'myseq.gtrig' not found
```

```r
# display
myseq.gtrig
```

```
## Error in eval(expr, envir, enclos): object 'myseq.gtrig' not found
```

Now we can perform a bootstrap, starting with our 'best' maximum likelihood tree, and performing nearest neighbour interchanges.


```r
myseq.gtrig.bs <- bootstrap.pml(myseq.gtrig,bs=100,trees=TRUE,optNni=TRUE)
```

```
## Error in bootstrap.pml(myseq.gtrig, bs = 100, trees = TRUE, optNni = TRUE): object 'myseq.gtrig' not found
```

As with the neighbour-joining tree, we can overlay bootstrap supports on our maximum likelihood tree.


```r
plotBS(myseq.gtrig$tree,myseq.gtrig.bs,type="phylogram",cex=0.5)
```

```
## Error in inherits(phy, "phylo"): object 'myseq.gtrig' not found
```

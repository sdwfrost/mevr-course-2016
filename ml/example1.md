# Model choice exercise



Load the libraries.


```r
library(ape)
library(phangorn)
```

Load in the sequence data and the neighbour joining tree.


```r
myalignment <- read.dna("ray2000_aligned.fas",format="fasta",as.matrix=TRUE)
mytree  <- nj(dist.dna(myalignment,"TN93"))
```

Now compare different models.


```r
myalignment.phydat <- as.phyDat(myalignment) # Convert to a format modeltest understand
myalignment.modeltest <- modelTest(myalignment.phydat,
                             tree=mytree,
                             model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"),
                             G = TRUE,
                             I = TRUE,
                             k = 4,
                             control = pml.control(epsilon = 1e-08, maxit = 3, trace = 1),
                             multicore = FALSE)
```

```
## Warning in pml(tree, data): negative edges length changed to 0!
```

```
## [1] "JC+I"
## [1] "JC+G"
## [1] "JC+G+I"
## [1] "F81+I"
## [1] "F81+G"
## [1] "F81+G+I"
## [1] "K80+I"
## [1] "K80+G"
## [1] "K80+G+I"
## [1] "HKY+I"
## [1] "HKY+G"
## [1] "HKY+G+I"
## [1] "SYM+I"
## [1] "SYM+G"
## [1] "SYM+G+I"
## [1] "GTR+I"
## [1] "GTR+G"
## [1] "GTR+G+I"
```

```r
myalignment.modeltest
```

```
##      Model  df    logLik      AIC      BIC
## 1       JC 139 -9093.606 18465.21 19023.80
## 2     JC+I 140 -8618.051 17516.10 18078.71
## 3     JC+G 140 -8368.362 17016.72 17579.33
## 4   JC+G+I 141 -8351.732 16985.46 17552.09
## 5      F81 142 -9005.856 18295.71 18866.35
## 6    F81+I 143 -8525.408 17336.82 17911.47
## 7    F81+G 143 -8251.858 16789.72 17364.38
## 8  F81+G+I 144 -8238.720 16765.44 17344.12
## 9      K80 140 -8482.996 17245.99 17808.60
## 10   K80+I 141 -7980.468 16242.94 16809.56
## 11   K80+G 141 -7685.580 15653.16 16219.78
## 12 K80+G+I 142 -7667.595 15619.19 16189.83
## 13     HKY 143 -8440.617 17167.23 17741.89
## 14   HKY+I 144 -7916.114 16120.23 16698.91
## 15   HKY+G 144 -7636.388 15560.78 16139.45
## 16 HKY+G+I 145 -7616.750 15523.50 16106.20
## 17     SYM 144 -8445.450 17178.90 17757.58
## 18   SYM+I 145 -7943.376 16176.75 16759.45
## 19   SYM+G 145 -7649.702 15589.40 16172.10
## 20 SYM+G+I 146 -7631.733 15555.47 16142.18
## 21     GTR 147 -8436.124 17166.25 17756.98
## 22   GTR+I 148 -7903.426 16102.85 16697.60
## 23   GTR+G 148 -7622.485 15540.97 16135.72
## 24 GTR+G+I 149 -7601.719 15501.44 16100.21
```

Using R, we can look for the model e.g. with the lowest AIC.


```r
myalignment.modeltest$Model[myalignment.modeltest$AIC==min(myalignment.modeltest$AIC)]
```

```
## [1] "GTR+G+I"
```

Now obtain a maximum likelihood tree, starting from the neighbour joining tree, and using the model chosen from ```modelTest```.


```r
myalignment.pml <- pml(mytree,myalignment.phydat,model="GTR+I+G",k=4)
```

```
## Warning in pml(mytree, myalignment.phydat, model = "GTR+I+G", k = 4):
## negative edges length changed to 0!
```

```r
# optimise
myalignment.pml <- optim.pml(myalignment.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)
```

```
## optimize edge weights:  -8617.809 --> -8415.953 
## optimize base frequencies:  -8415.953 --> -8306.962 
## optimize rate matrix:  -8306.962 --> -7739.162 
## optimize invariant sites:  -7739.162 --> -7662.088 
## optimize shape parameter:  -7662.088 --> -7662.088 
## optimize edge weights:  -7662.088 --> -7637.825 
## optimize topology:  -7637.825 --> -7621.82 
## optimize topology:  -7621.82 --> -7617.04 
## optimize topology:  -7617.04 --> -7606.44 
## 20 
## optimize base frequencies:  -7606.44 --> -7590.913 
## optimize rate matrix:  -7590.913 --> -7579.069 
## optimize invariant sites:  -7579.069 --> -7579.068 
## optimize shape parameter:  -7579.068 --> -7578.92 
## optimize edge weights:  -7578.92 --> -7577.282 
## optimize topology:  -7577.282 --> -7565.306 
## optimize topology:  -7565.306 --> -7557.973 
## optimize topology:  -7557.973 --> -7554.848 
## 15 
## optimize base frequencies:  -7554.848 --> -7551.863 
## optimize rate matrix:  -7551.863 --> -7550.22 
## optimize invariant sites:  -7550.22 --> -7550.168 
## optimize shape parameter:  -7550.168 --> -7550.117 
## optimize edge weights:  -7550.117 --> -7549.979 
## optimize topology:  -7549.979 --> -7542.122 
## optimize topology:  -7542.122 --> -7541.302 
## optimize topology:  -7541.302 --> -7541.302 
## 6 
## optimize base frequencies:  -7541.302 --> -7540.721 
## optimize rate matrix:  -7540.721 --> -7540.306 
## optimize invariant sites:  -7540.306 --> -7540.295 
## optimize shape parameter:  -7540.295 --> -7540.288 
## optimize edge weights:  -7540.288 --> -7540.272 
## optimize topology:  -7540.272 --> -7540.272 
## 0 
## optimize base frequencies:  -7540.272 --> -7540.091 
## optimize rate matrix:  -7540.091 --> -7539.995 
## optimize invariant sites:  -7539.995 --> -7539.994 
## optimize shape parameter:  -7539.994 --> -7539.992 
## optimize edge weights:  -7539.992 --> -7539.991 
## optimize base frequencies:  -7539.991 --> -7539.936 
## optimize rate matrix:  -7539.936 --> -7539.905 
## optimize invariant sites:  -7539.905 --> -7539.904 
## optimize shape parameter:  -7539.904 --> -7539.904 
## optimize edge weights:  -7539.904 --> -7539.904 
## optimize base frequencies:  -7539.904 --> -7539.887 
## optimize rate matrix:  -7539.887 --> -7539.878 
## optimize invariant sites:  -7539.878 --> -7539.878 
## optimize shape parameter:  -7539.878 --> -7539.878 
## optimize edge weights:  -7539.878 --> -7539.877 
## optimize base frequencies:  -7539.877 --> -7539.872 
## optimize rate matrix:  -7539.872 --> -7539.869 
## optimize invariant sites:  -7539.869 --> -7539.869 
## optimize shape parameter:  -7539.869 --> -7539.869 
## optimize edge weights:  -7539.869 --> -7539.869 
## optimize base frequencies:  -7539.869 --> -7539.868 
## optimize rate matrix:  -7539.868 --> -7539.867 
## optimize invariant sites:  -7539.867 --> -7539.867 
## optimize shape parameter:  -7539.867 --> -7539.867 
## optimize edge weights:  -7539.867 --> -7539.867 
## optimize base frequencies:  -7539.867 --> -7539.867 
## optimize rate matrix:  -7539.867 --> -7539.866 
## optimize invariant sites:  -7539.866 --> -7539.866 
## optimize shape parameter:  -7539.866 --> -7539.866 
## optimize edge weights:  -7539.866 --> -7539.866 
## optimize base frequencies:  -7539.866 --> -7539.866 
## optimize rate matrix:  -7539.866 --> -7539.866 
## optimize invariant sites:  -7539.866 --> -7539.866 
## optimize shape parameter:  -7539.866 --> -7539.866 
## optimize edge weights:  -7539.866 --> -7539.866 
## optimize base frequencies:  -7539.866 --> -7539.866 
## optimize rate matrix:  -7539.866 --> -7539.866 
## optimize invariant sites:  -7539.866 --> -7539.866 
## optimize shape parameter:  -7539.866 --> -7539.866 
## optimize edge weights:  -7539.866 --> -7539.866
```

```r
# display
myalignment.pml
```

```
## 
##  loglikelihood: -7539.866 
## 
## unconstrained loglikelihood: -1990.566 
## Proportion of invariant sites: 0.270931 
## Discrete gamma model
## Number of rate categories: 4 
## Shape parameter: 0.905875 
## 
## Rate matrix:
##          a          c         g         t
## a 0.000000  1.3729288 8.0371277  1.684299
## c 1.372929  0.0000000 0.7245706 11.250734
## g 8.037128  0.7245706 0.0000000  1.000000
## t 1.684299 11.2507341 1.0000000  0.000000
## 
## Base frequencies:  
## 0.1764514 0.3412765 0.2521093 0.2301627
```

This shows the maximum likelihood parameter values. The ML tree is contained in the fit as well.


```r
myalignment.pml$tree
```

```
## 
## Phylogenetic tree with 71 tips and 69 internal nodes.
## 
## Tip labels:
## 	AF271819, AF271820, AF271821, AF271823, AF271824, AF271822, ...
## 
## Unrooted; includes branch lengths.
```

We can compare our original tree and our new tree numerically using ```treedist```.


```r
treedist(myalignment.pml$tree,mytree)
```

```
##      symmetric.difference   branch.score.difference 
##                 56.000000                  1.024206 
##           path.difference quadratic.path.difference 
##                153.701008                 42.384989
```

This isn't particularly helpful. A comparison of edge lengths is quite striking.


```r
sum(myalignment.pml$tree$edge.length)
```

```
## [1] 7.13001
```

```r
sum(mytree$edge.length)
```

```
## [1] 3.115848
```

Now we can perform a bootstrap, starting with our 'best' maximum likelihood tree, and performing nearest neighbour interchanges.


```r
myalignment.pml.bs <- bootstrap.pml(myalignment.pml,bs=100,trees=TRUE,optNni=TRUE)
```

As with the neighbour-joining tree, we can overlay bootstrap supports on our maximum likelihood tree.


```r
plotBS(myalignment.pml$tree,myalignment.pml.bs,type="phylogram",cex=0.5)
```

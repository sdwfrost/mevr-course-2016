# Alignment example

## Prerequisites

First, load the libraries you need. We will use ```Biostrings```, which will allow us to translate sequences, ```msa```, which allows us to align using several external packages (bundled into R for easy installation), and ```seqinr``` for more sequence manipulation.


```r
library(Biostrings)
library(msa)
library(seqinr)
library(magrittr)
source("utils.R")
```

It's also useful to set any run-specific filenames at the beginning, so you can reuse the document more easily. For this example, I'll use the Ray et al. (2000) dataset, edited so that it isn't aligned.


```r
nucseq.fn <- "ray2000_edited.fas"
myref.fn <- "ray2000_ref.fas"
```

I'll use additional file stubs to keep track of what we've done to the original data.

## Align in nucleotide space

First, we will align in nucleotide space. This is advisable for non-coding data (e.g. the LTR of many viruses), and can be sufficient for many coding regions.


```r
nucseq <- readDNAStringSet(nucseq.fn)
nucseq
```

```
##   A DNAStringSet instance of length 71
##      width seq                                         names               
##  [1]   410 TCTTGGCTCTGTTGTCCTGT...TCTATCCTGGCCACGTAACT AF271819
##  [2]   410 TTTTGGCACTTCTCTCATGC...TTTATGCAGGCCATGTTACC AF271820
##  [3]   410 TTCTGGCCCTTCTTTCATGC...TCTATGCAGGCCATATCACC AF271821
##  [4]   410 TTCTGGCCCTTCTCTCATGC...TCTATCCAGGCCACATTACC AF271823
##  [5]   410 TTCTGGCCCTTCTCTCATGC...TCTATTCAGGCCACGTTACG AF271824
##  ...   ... ...
## [67]   411 CTCTTGGCACTTCTCTCGTG...TCTATACAGGGCACATCACC AF271868
## [68]   411 CTCTTGGCACTTCTCTCGTG...TCTATACGGGGCACATTACT AF271826
## [69]   411 CTCTTGGCACTTCTCTCTTG...TCTATACAGGGCACATCACT AF271837
## [70]   411 CTCTTGGCACTTCTCTCGTG...TCTACACAGGGCACATCACT AF271846
## [71]   411 CTTTTGGCACTTCTCTCGTG...TCTATACGGGGCACATCACT AF271840
```


```r
nucseq.align <- msa(nucseq,method="Muscle")
```

```
## Error in is(inputSeq, "character"): object 'nucseq' not found
```

Save the resulting alignment to a file.


```r
writeXStringSet(nucseq.align@unmasked,paste(nucseq.fn,".align",sep=""))
```

```
## Error in is(x, "XStringSet"): object 'nucseq.align' not found
```

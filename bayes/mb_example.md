# Getting a timetree with MrBayes

## Load libraries


```r
library(ape)
library(stringr)
```

## Set filenames


```r
myseqs.fn <- "H5N1.fas" # input file
myseqs.nex.fn <- "H5N1.nex" # output in Nexus format
myseqs.mb <- "H5N1_mb.nex" # output in Nexus with MrBayes block
```

## Read in data

The test data are of influenza A H5N1 sequences from birds sampled in China.


```r
myseqs <- read.dna(myseqs.fn,format="fasta",as.matrix=FALSE)
```

## Preprocess

Get sequence names and tip dates.


```r
myseqs.names <- names(myseqs)
myseqs.names
```

```
##  [1] "A_chicken_Fujian_1042_2005"    "A_duck_Fujian_897_2005"       
##  [3] "A_duck_Fujian_13_2002"         "A_swine_Fujian_F1_2001"       
##  [5] "A_Chicken_Guangdong_810_2005"  "A_goose_Guangdong_2216_2005"  
##  [7] "A_chicken_Guangdong_174_2004"  "A_chicken_Guangdong_191_2004" 
##  [9] "A_duck_Guangdong_4610_2003"    "A_duck_Guangdong_01_2001"     
## [11] "A_duck_Guangdong_12_2000"      "A_Goose_Guangdong_1_1996"     
## [13] "A_duck_Guangxi_793_2005"       "A_goose_Guangxi_345_2005"     
## [15] "A_chicken_Guangxi_2439_2004"   "A_chicken_Guangxi_2461_2004"  
## [17] "A_duck_Guangxi_1378_2004"      "A_duck_Guangxi_1681_2004"     
## [19] "A_duck_Guangxi_2291_2004"      "A_duck_Guangxi_351_2004"      
## [21] "A_goose_Guangxi_1097_2004"     "A_goose_Guangxi_2112_2004"    
## [23] "A_goose_Guangxi_914_2004"      "A_duck_Guangxi_22_2001"       
## [25] "A_duck_Guangxi_50_2001"        "A_greyheron_HongKong_837_2004"
## [27] "A_bird_HongKong_213_2003"      "A_chicken_HongKong_863_2002"  
## [29] "A_bird_HongKong_97_1998"       "A_chicken_HongKong_258_1997"  
## [31] "A_chicken_HongKong_915_1997"   "A_Goose_HongKong_w355_1997"   
## [33] "A_bird_HongKong_481_1997"      "A_bird_HongKong_485_1997"     
## [35] "A_bird_HongKong_488_1997"      "A_bird_HongKong_503_1997"     
## [37] "A_bird_HongKong_514_1997"      "A_bird_HongKong_532_1997"     
## [39] "A_bird_HongKong_542_1997"      "A_duck_Hunan_1265_2005"       
## [41] "A_duck_Hunan_139_2005"         "A_duck_Hunan_157_2005"        
## [43] "A_duck_Hunan_182_2005"
```

In this case, the tip dates are the last four letters in the sequence name.


```r
myseqs.tipdates <- as.integer(str_sub(myseqs.names,-4))
myseqs.tipdates
```

```
##  [1] 2005 2005 2002 2001 2005 2005 2004 2004 2003 2001 2000 1996 2005 2005
## [15] 2004 2004 2004 2004 2004 2004 2004 2004 2004 2001 2001 2004 2003 2002
## [29] 1998 1997 1997 1997 1997 1997 1997 1997 1997 1997 1997 2005 2005 2005
## [43] 2005
```

MrBayes requires Nexus format, with an added block giving instructions to MrBayes. We first save the data as Nexus format, and read back in to manipulate further.


```r
write.nexus.data(as.character(myseqs),myseqs.nex.fn,interleaved=FALSE,gap="-",missing="N")
myseqs.nex <- readLines(myseqs.nex.fn)
```

We then write part of the MrBayes block that specifies the model.


```r
mbblock1 <- "
begin mrbayes;
  set autoclose=yes;
  lset nst=6 rates=invgamma ngammacat=4;
	prset statefreqpr=fixed(empirical) brlenspr=clock:uniform;
  prset clockratepr=normal(0.00003,0.00001);
  prset treeagepr=exponential(0.0001);
  prset nodeagepr=calibrated;
"
```

We calculate the ages of the tips (rather than their times), and generate another part of the block that includes the calibration points. Note that we can have ranges on the sampling times (in this case, as we have sampling to the nearest year).


```r
myseqs.tipages <- max(myseqs.tipdates)-myseqs.tipdates
numseqs <- length(myseqs.tipages)
mbblock2 <- rep("",numseqs)
for(i in 1:numseqs){
    tipname <- myseqs.names[i]
    tipage <- myseqs.tipages[i]
    mbblock2[i] <- paste("calibrate ",tipname,"=Uniform(",tipage,",",tipage+1,");",sep="")
}
mbblock2 <- paste(mbblock2,collapse="\n")
```

We set the MCMC parameters.


```r
mbblock3 <- "
  mcmc ngen=1000000 nruns=2 nchains=2 samplefreq=1000;
end;
"
```

We then paste the blocks together and write to a file.


```r
myseqs.nexus.withmb <- paste(paste(myseqs.nex,collapse="\n"),mbblock1,mbblock2,mbblock3,sep="")
write(myseqs.nexus.withmb,file=myseqs.mb)
```

We can now run MrBayes with the command ```mb``` with the output filename as the single argument.

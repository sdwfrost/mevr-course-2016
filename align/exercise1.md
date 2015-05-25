# Alignment exercise

The aim of this objective is to illustrate the potential difficulties in aligning highly variable sequences, such as those from RNA viruses. I have provided two simulated sequence datasets generated using the program INDELible (Fletcher and Yang, 2009), for which the true alignment is known.

1. ```denv.fas```: simulated dengue virus sequences
2. ```influenza_HA.fas```: simulated influenza A haemagglutinin sequences

Align these as best you can, and then we will compare your answers to the true alignment. For your convenience, I have also provided translated sequences.

You may use the following commands in R (or any other alignment program you like). Using e.g. ```?msa``` to get more information on each command.

- ```msa``` from the library ```msa```
- ```mafft``` from the library ```ips```
- ```prank``` from the library ```ips```
- ```clustal``` (using either ```clustalw``` or Clustal-Omega) from the library ```ape```
- ```t_coffee``` from the library ```ape```

Examples of how to use these commands for nucleotide sequences are given below.

## Nucleotide alignment

### ```mafft``` and ```prank```

MAFFT is one of the fastest multiple sequence alignment programs available, and can be used for thousands of sequences.


```r
library(ips)
myseqs <- read.fas("denv.fas")
myseqs.align.mafft <- mafft(myseqs,path="/usr/bin/mafft")
write.fas(myseqs.align.mafft,"denv.fas.mafft")
```

Be warned, PRANK can take a long time...


```r
prank(myseqs,"denv.fas.prank",path="/usr/local/bin/prank")
```

### ```clustalw``` and ```clustalo```


```r
library(ape)
myseqs.clustalw <- clustal(myseqs,exec="/usr/bin/clustalw")
write.fas(myseqs.clustalw,"denv.fas.clustalw")
```

The ```clustal``` command in ```ape`` will also work for Clustal Omega (```clustalo```).


```r
myseqs.clustalo <- clustal(myseqs,exec="/usr/bin/clustalo")
write.fas(myseqs.clustalo,"denv.fas.clustalo")
```


```r
myseqs.tcoffee <- tcoffee(myseqs)
write.fas(myseqs.tcoffee,"denv.fas.tcoffee")
```

## Amino acid alignment

Try aligning in amino acid space too. Unfortunately, R does not have good support for amino acid sequences, so this is best done using system calls (which you can obviously also use for the nucleotide sequences too).

### ```mafft```


```r
cmd <- "mafft denv.fas.aa > denv.fas.aa.mafft"
system(cmd)
```

### ```muscle```


```r
cmd <- "muscle -in denv.fas.aa -out denv.fas.aa.muscle"
system(cmd)
```

### Clustal Omega


```r
cmd <- "clustalo -i denv.fas.aa -o denv.fas.aa.clustalo"
system(cmd)
```

### PRANK


```r
cmd <- "prank -d=denv.fas.aa -o=denv.fas.aa.prank"
system(cmd)
```


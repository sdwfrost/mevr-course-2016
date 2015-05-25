# Practice

There are a number of alignment packages that can be called from R.

- ```msa```: calls ClustalW, ClustalOmega, and MUSCLE (BioConductor package)
- ```ips```: calls mafft and prank.
- ```ape```: calls ClustalW, MUSCLE, T-Coffee

In addition, alignment programs can be called via ```system``` calls. I will go through a worked example, showing how difficult it can sometimes be to align highly variable sequences, followed by an exercise involving simulated data, for which the true alignments are known.

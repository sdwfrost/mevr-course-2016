# Retrieving sequences example






```r
library(XML)
library(rentrez)
library(ape)
library(phangorn)
```

```
## Warning: closing unused connection 5
## (ray2000_edited.fas.mapped.degapped.aa)
```

```r
source("utils.R")
```

You also need to set your email address.


```r
entrez_email <- "sdwfrost@gmail.com"
```

# Getting sequences e.g. from a paper

## Basic

If we have a list of accessions, we can post these to the nucleotide database, known as ```nuccore```, and then get the sequences back in one of several formats.

You can find the list of search terms using ```entrez_db_searchable```.


```r
search_fields <- entrez_db_searchable("nuccore")
search_fields
```

```
## Searchable fields for database 'nuccore'
##  [1] ALL  UID  FILT WORD TITL KYWD AUTH JOUR VOL  ISS  PAGE ORGN ACCN PACC
## [15] GENE PROT ECNO PDAT MDAT SUBS PROP SQID GPRJ SLEN FKEY PORG COMP ASSM
## [29] DIV  STRN ISOL CULT BRD  BIOS
```



```r
myaccs <- c("AF271819", "AF271820", "AF271821", "AF271823", "AF271824", 
"AF271822", "AF271818", "AF271817", "AF271876", "AF271875", "AF271882", 
"AF271879", "AF271878", "AF271877", "AF271887", "AF271881", "AF271880", 
"AF271874", "AF271883", "AF271885", "AF271886", "AF271884", "AF271847", 
"AF271851", "AF271866", "AF271839", "AF271841", "AF271849", "AF271873", 
"AF271843", "AF271858", "AF271867", "AF271855", "AF271872", "AF271831", 
"AF271825", "AF271863", "AF271835", "AF271845", "AF271859", "AF271860", 
"AF271836", "AF271832", "AF271827", "AF271834", "AF271854", "AF271828", 
"AF271850", "AF271853", "AF271848", "AF271857", "AF271829", "AF271865", 
"AF271838", "AF271856", "AF271852", "AF271864", "AF271862", "AF271842", 
"AF271830", "AF271871", "AF271869", "AF271833", "AF271844", "AF271870", 
"AF271861", "AF271868", "AF271826", "AF271837", "AF271846", "AF271840"
)
myids <- accToIds(myaccs)
mypost <- entrez_post(db="nuccore",id=myids)
myseqs <- entrez_fetch(db="nuccore", ids="", rettype="fasta", WebEnv=mypost$WebEnv, query_key=mypost$QueryKey)
cat(myseqs,file="query.fas")
```

## Advanced

We search PubMed for the paper identifier.


```r
ray2000 <- entrez_search(db="pubmed",term="ray[au] and hcv and egypt and 2000", retmax=10,usehistory=TRUE)
ray2000
```

```
## Entrez search result with 1 hits (object contains 1 IDs and a cookie)
```

We can now get the data in PubMed on this - very useful for literature reviews!


```r
ray2000.pubmed <- entrez_fetch(db="pubmed",id=ray2000$ids,rettype="xml")
```


```r
ray2000.link <- entrez_link(db="nuccore",dbfrom="pubmed",linkname="pubmed_nuccore",id=ray2000$ids)
```


```r
ray2000.ids <- ray2000.link$pubmed_nuccore
ray2000.post <- entrez_post(db="nuccore",id=ray2000.ids)
```

```
## Error in function (type, msg, asError = TRUE) : Recv failure: Connection reset by peer
```

```r
ray2000.seqs <- entrez_fetch(db="nuccore", ids="", rettype="fasta", WebEnv=ray2000.post$WebEnv, query_key=ray2000.post$QueryKey)
```

```
## Error in make_entrez_query("efetch", db = db, rettype = rettype, config = config, : object 'ray2000.post' not found
```

```r
ray2000.xml <- entrez_fetch(db="nuccore", ids="", rettype="xml", WebEnv=ray2000.post$WebEnv, query_key=ray2000.post$QueryKey)
```

```
## Error in make_entrez_query("efetch", db = db, rettype = rettype, config = config, : object 'ray2000.post' not found
```


```r
write(ray2000.seqs,"ray2000.fas")
```

```
## Error in cat(x, file = file, sep = c(rep.int(sep, ncolumns - 1), "\n"), : object 'ray2000.seqs' not found
```

# Getting all genomes from a virus species

To get all genomes manually we:

- Search NCBI Taxonomy for the species and retrieve the taxonomic identifier
- Link to the genome database to get the reference genome ID
- Link to the nucleotide database to get the nucleotide identifiers for other genomes of the species
- Post the nucleotide IDs to the database
- Retrieve sequences in a specific format

## Searching the taxonomy database with some text


```r
tax <- entrez_search(db="taxonomy",term="Hepatitis C", retmax=10,usehistory=TRUE)
tax
```

```
## Entrez search result with 1 hits (object contains 1 IDs and a cookie)
```

```r
tax$ids
```

```
## [1] "11103"
```

## Using the taxonomic ID to search for genomes


```r
genome <- entrez_search(db="genome",term=paste("txid",tax$ids,"[Orgn]",sep=""))
genome
```

```
## Entrez search result with 1 hits (object contains 1 IDs and no cookie)
```

There is a lot of information associated with the genome, that we can obtain using ```entrez_summary```.


```r
genome.summary <- entrez_summary(db="genome",id=genome$ids)
genome.summary
```

```
## esummary result with 21 items:
##  [1] uid                   organism_name         organism_kingdom     
##  [4] organism_group        organism_subgroup     defline              
##  [7] projectid             project_accession     status               
## [10] number_of_chromosomes number_of_plasmids    number_of_organelles 
## [13] assembly_name         assembly_accession    assemblyid           
## [16] create_date           options               weight               
## [19] chromosome_assemblies scaffold_assemblies   sra_genomes
```

## Linking from genome to nucleotide databases

To cross reference databases - like clicking a link in the table - we use ```entrez_link```.


```r
link <- entrez_link(db="nuccore",dbfrom="genome",linkname="genome_nuccore_samespecies",id=genome$ids)
link
```

```
## elink result with ids from 1 databases:
## [1] genome_nuccore_samespecies
```

```r
seq.ids <- link$genome_nuccore_samespecies
length(seq.ids)
```

```
## [1] 1430
```

## Getting sequences given a list of sequence identifiers

Now that we have the sequence identifiers we need, we can download the sequences. I'll just get 100 genomes. Firstly, we post the IDs to the database using ```entrez_post```.


```r
seqs.post <- entrez_post(db="nuccore",id=seq.ids[1:100])
```

Next, we use ```entrez_fetch``` to fetch the sequences.


```r
seqs.fasta <- entrez_fetch(db="nuccore", ids="", rettype="fasta", WebEnv=seqs.post$WebEnv, query_key=seqs.post$QueryKey, retend=100)
```

```
## Error in function (type, msg, asError = TRUE) : Recv failure: Connection reset by peer
```


```r
write(seqs.fasta,"hcv_genomes.fas")
```

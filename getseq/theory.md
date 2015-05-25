# Theory

- The first thing you will need in an analysis are some data!
- These may come from a variety of sources, and in a variety of formats

## Public databases

- General databases, part of the International Nucleotide Sequence Database Collaboration
  - [GenBank](http://www.ncbi.nlm.nih.gov/genbank/)
  - [European Nucleotide Archive (ENA)](http://www.ebi.ac.uk/ena)
  - [DNA Databank of Japan (DDBJ)](http://www.ddbj.nig.ac.jp/)
- Virus specific databases, e.g.
  - [NCBI Influenza Virus Resource](http://www.ncbi.nlm.nih.gov/genomes/FLU/FLU.html)
  - [Los Alamos HIV Sequence Database](http://hiv-web.lanl.gov)

## Genbank

- We will download sequences and associated data from GenBank
- Central to this goal is the way in which sequences are identified
- GI numbers are simply a series of digits that are assigned consecutively
- Accession numbers in the nucleotide database comprise of 1-2 letters, 5-6 digits, a dot and a version number
- In the reference sequence database, accession numbers comprise of 2 letters, an underscore, and several digits

## Annotations

- At a minimum, GenBank expects the following for each viral sequence submitted
  - Organism name
  - Isolate
  - Country (optional)
  - Collection date (optional)
  - Host (optional)
- For influenza virus, rotavirus, and norovirus, the above optional fields are mandatory

## Entrez through the web

- Multiple databases at NCBI can be searched via the search engine Entrez
- For exploratory searches, it is useful to access Entrez directly through a web browser at [http://www.ncbi.nlm.nih.gov/Entrez](http://www.ncbi.nlm.nih.gov/Entrez)

## Output formats

- There are a number of different formats for output once sequences have been identified; we will consider three:
  - FASTA
  - List of GI numbers
  - INSDSeq XML

## Entrez through R

- Eutils can also be used to access Entrez
- This has the advantages of being able to reproduce the searches, as well as update them, as new data become available

## Extracting annotation data

- Information like host and collection data is not provided in the FASTA file; it has to be extracted from the full record
- In the practical, we will see how this information can be extracted from the record in INSDSeq XML format

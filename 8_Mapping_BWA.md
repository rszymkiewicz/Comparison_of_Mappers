# Mapping using BWA
## Commands and Parameters  
### ***Library Preparation***  
```
bwa index [reference file in fasta format]
```  
### ***Mapping***
Alignment using *bwa mem* uses seeding alignments based upon the maximal exact matches (MEM) and extends seeds using the affine-gap Smith-Waterman algorith (SW). The algorithm uses *local alignment* having the ability to produce multiple alignments for certain regions of the reference sequence. 
```
bwa mem -t [reference file in fasta format] [Forward paired strand in fastq format] [Reverse paired strand in fastq format] > [Output file in sam format]
``` 
Parameter | Description
----------|------------
-t | *flag indicating the number of threads*  

**BWA-MEM manual:** http://bio-bwa.sourceforge.net/bwa.shtml

## Python Script

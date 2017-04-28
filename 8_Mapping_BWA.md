# Mapping using BWA
**BWA Manual:** http://bio-bwa.sourceforge.net/bwa.shtml  
Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]  
Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].  

BWA stands for Burrows-Wheeler Alignment Tool and has a total of 3 algorithms to which can be used to align a read to a reference sequence, BWA-backtrack, BWA-SW, and BWA-MEM.  

**BWA-MEM** is suggested for Illumina reads which range in size between 70 and 100bp and provides both fast and accurate high-quality results. BWA-MEM has the ability to tolerate a higher error rate due to the capability to predict longer alignments. Error rates include 2% for a 100bp alignment, 3% for 200bp alignment, 5% for 500bp alignment, and up to 10% for alignments of 1000bp or longer. As our sample set has Illumina sequences ranging in size we will use **BWA-MEM** algorithim to align our Illumina reads to our reference sequence. A python script used to perform BWA-MEM alignment of our sample set can be found at the bottom of the page. 

## Commands and Parameters  
### *Library Preparation*  
BWA requires an FM-index, otherwise known as a library, prior to performing alignment. **BWA-index** command is used to prepare the index database sequence library.  

```
bwa index [reference file in fasta format]
```  

### *Mapping*

BWA outputs one alignment per read however; if the read is chimeric it may output multiple alignments. 
BWA reports mapping quality based on each individual read compared to reporting mapping quality per pair when using paired-end features. THis is due to the ability of one strand to have the ability to map to a tandom repeat whereas the other strand does not therefore accurate mapping cannot be determined in these case scenarios. 

Alignment using *bwa mem* uses seeding alignments based upon the maximal exact matches (MEM) and extends seeds using the affine-gap Smith-Waterman algorith (SW). The algorithm uses *local alignment* having the ability to produce multiple alignments for certain regions of the reference sequence. 

```
bwa mem -t [reference file in fasta format] [Forward paired strand in fastq format] [Reverse paired strand in fastq format] > [Output file in sam format]
``` 
Parameter | Description
----------|------------
-t | *flag indicating the number of threads*  

## BWA-MEM Analysis - Sample P750    
***BWA-MEM Sam Output***  

# Python Script - Mapping with BWA

## Now that we have completed alignment of our samples to our reference sequence, let us continue and assess the possible SNPs and indels which maybe present within our sample set in comparison to our reference sequence. To do this lets learn about [*variant calling*](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/9_Variant_Calling.md).   

#### Return to *Bowtie2 Mapping* by clicking [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/7_Mapping_Bowtie2.md).  

# Mapping using BWA
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

lignment using *bwa mem* uses seeding alignments based upon the maximal exact matches (MEM) and extends seeds using the affine-gap Smith-Waterman algorith (SW). The algorithm uses *local alignment* having the ability to produce multiple alignments for certain regions of the reference sequence. 

```
bwa mem -t [reference file in fasta format] [Forward paired strand in fastq format] [Reverse paired strand in fastq format] > [Output file in sam format]
``` 
Parameter | Description
----------|------------
-t | *flag indicating the number of threads*  

## BWA-MEM Analysis - Sample P750    
***BWA-MEM Sam Output***  
Unlike Bowtie2 **BWA-MEM** does not produce a simplistic summary upon completion of alignment. In order to obtain a summary of information we can use **```samtools flagstat```** to aquire some alignment information from the **BWA-MEM** tool. To do this we first must convert our *sam* file to a *bam* file as outlined in the [next tutorial page]().  Below is a general outline of the **samtools flagstat** command followed by the output for Sample P750.   

```
samtools flagstat [BWA_bam_output_file]
```  
```
1521272 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
1521272 + 0 mapped (100.00%:-nan%)
1521272 + 0 paired in sequencing
761222 + 0 read1
760050 + 0 read2
1507304 + 0 properly paired (99.08%:-nan%)
1520163 + 0 with itself and mate mapped
1109 + 0 singletons (0.07%:-nan%)
83 + 0 with mate mapped to a different chr
9 + 0 with mate mapped to a different chr (mapQ>=5)
```  

A total of 1,521,272 reads successfully passed the quality threshold and were analyzed. 761,222 reads were from the forward strand and 760,050 were from the reverse strand. A total of 1,507,304 reads were mapped in a proper pair (99.08%) to the reference sequence. A total of 1109 reads were not mapped however their mate did map to the reference sequence. 83 paired reads mapped to a different chromosome with a total of 9 paired reads mapping to a different chromosome with a quality greater than 5.  

# Manual and Citation  
**BWA Manual:** http://bio-bwa.sourceforge.net/bwa.shtml  
Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]  
Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].  

# Python Script - Mapping with BWA
## [BWA-MEM](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/BWA_MEM.py)

## Now that we have completed mapping and obtained output files in *sam* file format, let convert the *sam* file format into [*other useful file formats*](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/9_sam_to_bam_to_sorted.md).   

# Variant Calling using *samtools* and *bedtools*  
**samtools manual:** http://www.htslib.org/doc/samtools.html  
**bcftools manual:** http://www.htslib.org/doc/bcftools.html  
Li H A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627] 

## samtools and bcftools
***samtools*** is a bioinformatic suite tool which manipulates files in SAM format. A SAM file stands for *sequence alignment/map* and contains nucleotide sequence alignments. The samtools suite provides a variety of functionality including sorting, merging, indexing, and generating alignments per nucleotide position. In many instances, Unix pipes are utilized in the samtools commands in order to combine various suite features across the htslib group of tools. In this case we take the output provided from the samtools command and use it as input for the bcftools command. 

***bcftools*** is another bioinformatic suite tool and is similar to samtools however; it has the ability to manipulate vcf file formats. In our case, the samtools command we will be using will output a vcf file and the bcftools command will manipulate that vcf file and allow us to extract specific information for further downstream analysis.  

nFor our purposes, we will be using the mpileup feature within the samtools suite in combination with the bcftools call command in order to create a pileup of reads which aligned to a specific reference sequence which can be used for variant calling analysis. To do this we will be taking a sorted bam file from each of the two mapping processes we previously completed [Bowtie2](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/7_Mapping_Bowtie2.md) and [BWA](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/8_Mapping_BWA.md) as input for our *mpileup* command. We will continue to use our Sample P750 as an example for this part of the tutorial and a python script will be used for the remaining samples of our sample set. The python script can be found at the bottom of the page.   

## Commands and Parameters  
### ***samtools Library Index***
```
samtools faidx [reference file in fasta format]
 ``` 
### ***samtools mpileup and bcftools call***  
**samtools mpileup** extracts information from an input file in BAM format (binary version of SAM file) and combines all the read alignment information into a VCF/BCF file format as output.  
**bcftools call** command calls the variant sites estimated by the samtools mpileup command and attempts to distinguish genotypes.  
e
```
samtools mpileup -vf [reference file in fasta format] [sample file in sorted bam format] | bcftools call -vmO z -o [output file in vcf gzipped format]
```  
#Tool | Command | Parameter | Description  
----|---------|-----------|------------  
 samtools | mpileup | -v | *variant calling format (vcf).*  
| | | -f | *flag indicating reference file in fasta format.*
| | | -t DP | *flag indicating the output file will include per-sample read depth.*  Note: This can only be used in conjungtion with -v or -g parameters of the mpileup suite.   
| | | -t SP | *flag indicating the output file will include per-sample Phred-scaled strand bias P-value.*  Note: This can only be used in conjungtion with -v or -g parameters of the mpileup suite.  bcftools | call | -v | *flag indicating to output variant sites only.*
| | | -m | *multiallelic and rare-variant calling feature of bcftools.*  
| | | -O z | *flag indicating the output file type; the 'z' represents an uncompressed VCF file format.*   
| | |-o | *flag indicating output file path and file format.*  

## File Output
### *General information pertaining to VCF/BCF File Output*
**VCF/BCF** is a variant call format file or binary call format file which contains the genotype likelihoods. There are 3 components to the VCF file: metadata header, header lines, and data lines. The header consists of information explaining the contents of each of the fields of the file represented by an inital hashtag (#) followed by a description. Each of the data lines of the vcf output file contains a number of features in a 10 TAB-delimited column format. It is important to note that column 9 and 10 are not always included within the vcf output file as it is dependent on the parameters used within the command. A chart depicting the 10 columns of the data line section can be seen below:  

Column Number | Column Name | Description  
--------------|-------------|---------------  
1 | CHROM | *Chromosome name.*  
2 | POS | *1-based position. The position prior to an indel is represented for indels.*  
3 | ID | *Variant identifier.*  
4 | REF | *Reference sequence position involved with the variant. A SNP is represented by a single base.*   
5 | ALT | *Comma list of alternative sequences.*  
6 | QUAL | *Phred-scaled probability of all samples being homozygous reference.*   
7 | FILTER | *Semicolon list of filters to which the variant fails.*  
8 | INFO | *Semicolon list describing variant information.*  
9 | FORMAT | *Colon list of indivdiual genotypes.*  
10 | SAMPLE | *Genotype information from FORMAT column.*  

Of particular interest is the **read base** column which contains important information pertaining to variant calling. These items include information pertaining to matches, mismatches, indels, strands, mapping quality, read start and read end and can be depicted in the chart below: 

VCF/BCF Feature | Description  
---------------|-------------  
. | ***match to reference*** *forward strand*.  
, | ***match to reference*** *reverse strand*.  
```> or <``` | ***reference skip***.  
ACGTN | ***mismatch to reference*** *forward strand*.  
acgtn | ***mismatch to reference*** *reverse strand*.  
```\\+[0-9]+[ACGTNacgtn]+``` | ***insertion between reference positons***.  
```\\+**[0-9]**+[ACGTNacgtn]+``` | *length of insertion represented by the integer provided.*   
```\\+[0-9]+**[ACGTNacgtn]**+``` | *identified insertion sequence.*  
`*` | ***deletion from the reference***.  
^ | *indicates the ***start*** of a read.*  
$ | *indicates the ***end*** of a read.*  
^ **#** | *indicates the mapping quality in ASCII format.*  Note: must subtract 33 from this integer to obtain the true mapping quality.  
## samtools and bedtools Analysis - Sample P750


# Python Script

## Lets investigate how we can compare our results from our two alignment tools (Bowtie2 and BWA) [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/11_Comparison.md).

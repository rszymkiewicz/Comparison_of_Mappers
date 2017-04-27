# Variant Calling using *samtools* and *bedtools*  
**samtools manual:** http://www.htslib.org/doc/samtools.html  
Li H A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627] 

## samtools 
***samtools*** is a bioinformatic suite tool which manipulates files in SAM format. A SAM file stands for *sequence alignment/map* and contains nucleotide sequence alignments. The samtools suite provides a variety of functionality including sorting, merging, indexing, and generating alignments per nucleotide position. In many instances, Unix pipes are utilized in the samtools commands in order to combine various suite features across the htslib group of tools.  

## bedtools
***bedtools***  



For our purposes, we will be using the mpileup feature within the samtools suite followed by the bedtools stats command in order to create a pileup of reads which aligned to a specific reference sequence  and statistical information which can be used for variant calling analysis.  

## Commands and Parameters  
 
### ***samtools Library Index***
```
samtools faidx [reference file in fasta format]
``` 
### ***samtools mpileup***  
**samtools mpileup** extracts information from an input file in BAM format (binary version of SAM file) and combines all the read alignment information into a VCF/BCF file format as output.  

```
samtools mpileup -vf [reference file in fasta format] [sample file in sorted bam format] | bcftools call -vmO z -o [output file in vcf gzipped format]
```  
Tool | Command | Parameter | Description  
----|---------|-----------|------------  
samtools | mpileup | -v | *variant calling format (vcf).*  
| | | -f | *flag indicating reference file in fasta format.*
| | | -t DP | *flag indicating the output file will include per-sample read depth.*  Note: This can only be used in conjungtion with -v or -g parameters of the mpileup suite.   
| | | -t SP | *flag indicating the output file will include per-sample Phred-scaled strand bias P-value.*  Note: This can only be used in conjungtion with -v or -g parameters of the mpileup suite.  
bcftools | call | -v |
| | | -m |   
| | | -O |   
| | | z |   
| | |-o | *flag indicating output file path and file format.*  

### ***bedtools stats***  
```
bcftools stats --threads 11 -F [reference file in fasta format] -s [sample file in sam format] > [output file and path]
```  
Parameter | Description
----------|------------  
-threads | *number of threads available for use.*  
-F | *flag indicating reference file in fasta format.*  
-s | *flag indicating sample file in sam format.*  

### ***bcftools/plot-vcfstats***  
```
bcftools/plot-vcfstats -p [output folder name] [sample file from *bcftools stats* output file]
```  
Parameter | Description
----------|------------
-p |  

## File Output
### *VCF/BCF File Output*
**VCF/BCF** is a variant call format file or binary call format file which contains the genotype likelihoods. Each line of the vcf/bcf output file contains a number of features in column format including: genomic position, chromosome name, 1-based coordinate reference base, number of reads covering the site, read bases, base qualities and alignment mapping qualities. Of particular interest is the **read base** column which contains important information pertaining to variant calling. These items include information pertaining to matches, mismatches, indels, strands, mapping quality, read start and read end and can be depicted in the chart below: 

VCF/BCF Feature | Description  
---------------|-------------  
. | ***match to reference*** *forward strand*.  
, | ***match to reference*** *reverse strand*.  
> or < | ***reference skip*** .  
ACGTN | ***mismatch to reference*** *forward strand*.  
acgtn | ***mismatch to reference*** *reverse strand*.  
\\+[0-9]+[ACGTNacgtn]+ | ***insertion between reference positons***  
\\+**[0-9]**+[ACGTNacgtn]+ | *length of insertion represented by the integer provided.*  
\\+[0-9]+**[ACGTNacgtn]**+ | *identified insertion sequence.*
`*` | ***deletion from the reference***. 
^ | *indicates the ***start*** of a read.*  
$ | *indicates the ***end*** of a read.*  
^ **#** | *indicates the mapping quality in ASCII format.*  Note: must subtract 33 from this integer to obtain the true mapping quality. 

## samtools and bedtools Analysis - Sample P750


# Python Script


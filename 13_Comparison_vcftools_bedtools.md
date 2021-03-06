#### Return to *Data Comparison Overview* by clicking [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/12_Comparison.md). 

# Data Comparison using *vcftools* and *bedtools*  
A variety of tools are available to assist in downstream analysis of Illumina sequencing samples. For our purposes, we will be using 3 bioinformatic suite tools, ***samtools, vcftools, bedtools,*** and ***bcftools***. Each bioinformatic tool has numerous commands available to assist in the manipulation of file outputs supporting the abundant information that can be obtained through genotyping.  

Our comparison of mapper tools (Bowtie2 and BWA) will begin with samtools using the *tabix* command to create a tab-deliminated genome position file. We will use the *compare* command of vcftools to output a venn diagram of the unique and similar results obtained from these two mapping tools. The *stats* command within the bcftools suite will produce various statistical information in a graphical format. Finally bedtools intersect and subtract commands will be used as another mechanism to identify overlaps and unique features when comparing the two mapping tools as used in this tutorial. 

## Commands and Parameters  
### samtools - *tabix*
**tabix** creates a compressed and indexed TAB-delimited file from common position sorted files. These files are required as input for *vcf-compare.*  

```
tabix -p vcf -h [sample file in vcf format]
```  
Parameter | Description
----------|-------------
-p | *flag indicating input sample file format [gff, bed, sam, vcf] required for indexing*  
-h | *print header and meta lines.*


### vcftools - *vcf-compare*  
**vcf-compare** compares the positions of multiple vcf files. In our case we will be comparing the Bowtie2 vcf file output to the BWA vcf file output to see if the two mapper tools produced any differences in variant calling.
```
vcf-compare [sample file in vcf format] [another sample in vcf format] > [output path]
```  
##### File Output  
Output information is based upon the first sample input vcf file and outlines the number of positions present within this file and absent within the other sample input vcf files. In addition, the perl script computes nonreference discordance rates known as multiallelic sites and compares actual sequence. 

### bcftools - *stats*
**stats** produces a text file of statistical information which can then be plotted using the plot-vcfstats command.  
```
bcftools stats --threads 11 -F [reference file in fasta format] -s [sample file in sam format] > [output file and path]
```  
Parameter | Description
----------|------------  
-threads | *number of threads available for use.*  
-F | *flag indicating reference file in fasta format.*  
-s | *flag indicating sample file in sam format.*

### bcftools - *plot-vcfstats*
**plot-vcfstats** is a python script which takes the bcftools stats output and creates various plots in a graphical format. A pdf format of all plots is also created as output with this command.   
```
bcftools/plot-vcfstats -p [output folder name] [sample file from *bcftools stats* output file]
```  
Parameter | Description
----------|------------
-p | *flag indicating the output file prefix.*   

### bedtools - *intersect*  
**intersect** compares sample files in order to identify the genomic components which overlap with one another.  
```
bedtools intersect -a [sample file in vcf format] -b [another sample in vcf format] > [output file and path]
```  
Parameter | Description
----------|------------
-a | *flag indicating the 1st input sample file in vcf format to be compared. The 1st input sample file ('a') is compared to the remaining input files in order to identify overlaps.*  
-b | *flag indicating the 2nd input sample file in vcf format to be compared.*  

##### File Output  
The bedtools intersect tool outputs a vcf file of common genome features. To review the description of a VCF file we can return to our [previous tutorial page](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/11_Variant_Calling_Samtools_Bcftools.md).

### bedtools - *subtract*  
**subtract** compares sample files in order to identify the genomic components which are unique. The second input sample file ('b') is compared to the 1st input sample file ('a') and all regions which overlap are disregarded and not reported. Only regions which remain from the 1st input sample file ('a') are reported as output.  
```
bedtools subtract -a [sample file in vcf format] -b [another sample file in vcf format] > [output file and path]
```  
Parameter | Description
----------|------------
-a | *flag indicating the 1st input sample file in vcf format to be compared.*  
-b | *flag indicating the 2nd input sample file in vcf format to be compared. The 2nd input sample file ('b') is used to identify unique features of the 1st input sample file ('a').*  

##### File Output 
The bedtools intersect tool outputs a vcf file of unique genome features based on the base file provided in the -a parameter. To review the description of a VCF file we can return to our [previous tutorial page](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/11_Variant_Calling_Samtools_Bcftools.md).

# Analysis of Sample P750
### Output files for this page of the tutorial can be found in the following links:
[vcf-compare](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/P750.compared)  
[bedtools intersect](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/Bedtools_Intersect_COMMON_B2_BWA.vcf)  
[bedtools subtract (Bowtie2 unique)](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/Unique_B2.vcf)  
[bedtools subtract (BWA unique)](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/Unique_BWA.vcf)  

***vcf-compare***  
```
VN	14	Sample_P750_BWA.vcf.gz (15.4%)
VN	28	Sample_P750_B2.vcf.gz (26.7%)
VN	77	Sample_P750_B2.vcf.gz (73.3%)	Sample_P750_BWA.vcf.gz (84.6%)
```  
Based on the ***vcf-compare*** output it can be determined that a total of 14 variant sites were unique to BWA mapping. A total of 28 variant sites were unique to Bowtie2 mapping. In addition, a total of 77 variant sites were shared between Bowtie2 and BWA mapping for Sample P750 which represented 73.3% and 84.6% of the total variant sites per mapping tool respectively.  

***bcftools stats and bcftools/plots-vcfstats***  
**Bowtie 2**  
- 84 SNPs identified
- 21 indels identified
- Ratio of transitions and transversions (ts/tv) 2.111111

**BWA**
- 74 SNPs identified
- 17 indels identified
- Ratio of transitions and transversions (ts/tv) 4.285714

***bcftools intersect*** 
A total of 78 regions were common between the two files.

***bcftools subtract*** 
A total of 14 regions were unique to BWA mapping consisting of 8 SNPs and 6 indels.  
A total of 27 regions were unique to Bowtie2 mapping consisting of 18 SNPs and 9 indels.   

## Manuals and Citation
**Tabix Manual from htslib-1.31:** http://www.htslib.org/doc/tabix.html  
Heng Li; Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 2011; 27 (5): 718-719. doi: 10.1093/bioinformatics/btq671  

**vcf-compare perl module:** https://vcftools.github.io/perl_module.html#vcf-compare  
The Variant Call Format and VCFtools, Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert Handsaker, Gerton Lunter, Gabor Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin and 1000 Genomes Project Analysis Group, Bioinformatics, 2011 http://dx.doi.org/10.1093/bioinformatics/btr330  

**bedtools intersect manual:** http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html  
**bedtools subtract manual:** http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html  
Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842. 

# Python Script - Variant Calling Comparison of Bowtie2 and BWA-MEM
## [Compare Variant Calling](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/CompareVariantCalling.py)  

### Important Note:  
Prior to using the *Compare_Variant_Calling.py* python script, one will need to run the following to scripts.  
[VCF File Paths (Bowtie2)](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/VCF_Path_Location_Bowtie2.py)  
[VCF File Paths (BWA-MEM)](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/VCF_Path_Location_BWA.py)  

Upon running the two scripts above, two text output files will be produced (one per script) representing the vcf output file paths from the variant calling tutorial page from each of the two mappers (Bowtie2 and BWA_MEM). The location of these two text files must be in the same location where the *Compare_Variant_Calling.py* python script is run.  

## Please click [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/14_Conclusion.md) for a summary of this tutorial.

 # Data Comparison using *vcftools* and *bedtools*  

## Commands and Parameters  
### Tabix 
```
tabix -p vcf [sample file in vcf format]
```  
Parameter | Description
----------|-------------
-p | *flag indicating input sample file format [gff, bed, sam, vcf] required for indexing*  

**Tabix Manual from htslib-1.31:** http://www.htslib.org/doc/tabix.html

### Vcftools  
#### ***vcf-compare***  
Compares the positions of multiple vcf files. 
```
vcf-compare [sample file in vcf format] [another sample in vcf format] > [output path]
```  
##### File Output  
Output information is based upon the first sample input vcf file and outlines the number of positions present within this file and absent within the other sample input vcf files. In addition, the perl script computes nonreference discordance rates known as multiallelic sites and compares actual sequence. 

**vcf-compare perl module:** https://vcftools.github.io/perl_module.html#vcf-compare  

### Bedtools  
#### ***intersect***
Compares sample files in order to identify the genomic components which overlap with one another.  
```
bedtools intersect -a [sample file in vcf format] -b [another sample in vcf format] > [output file and path]
```  
Parameter | Description
----------|------------
-a | *flag indicating the 1st input sample file in vcf format to be compared. The 1st input sample file ('a') is compared to the remaining input files in order to identify overlaps.*  
-b | *flag indicating the 2nd input sample file in vcf format to be compared.*  

##### File Output  
The bedtools intersect tool outputs a vcf file of common genome features.

**bedtools manual:** http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html  

#### ***subtract***  
Compares sample files in order to identify the genomic components which are unique. The second input sample file ('b') is compared to the 1st input sample file ('a') and all regions which overlap are disregarded and not reported. Only regions which remain from the 1st input sample file ('a') are reported as output.  
```
bedtools subtract -a [sample file in vcf format] -b [another sample file in vcf format] > [output file and path]
```  
Parameter | Description
----------|------------
-a | *flag indicating the 1st input sample file in vcf format to be compared.*  
-b | *flag indicating the 2nd input sample file in vcf format to be compared. The 2nd input sample file ('b') is used to identify unique features of the 1st input sample file ('a').*  
##### File Output 
The bedtools intersect tool outputs a vcf file of unique genome features based on the base file provided in the -a parameter.

**bedtools manual:** http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html

## Python Script  

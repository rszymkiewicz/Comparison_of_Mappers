# Data Comparison using *vcftools* and *bedtools*  

## Commands and Parameters  
### Tabix 
```
tabix -p vcf [sample file in vcf format]
```  
Parameter | Description
-p | **  

### Vcftools  
#### ***vcf-compare***  
```
vcf-compare [sample file in vcf format] [another sample in vcf format] > [output path]
```  
##### File Output  
  

### Bedtools  
#### ***intersect***  
```
bedtools intersect -a [sample file in vcf format] -b [another sample in vcf format] > [output file and path]
```  
Parameter | Description
-a | *flag indicating the 1st sample file in vcf format to be compared*  
-b | *flag indicating the 2nd sample file in vcf format to be compared*  

##### File Output  
The bedtools intersect tool outputs a vcf file of common genome features.
#### ***subtract***  
```
bedtools subtract -a [sample file in vcf format] -b [another sample file in vcf format] > [output file and path]
```  
Parameter | Description
-a | *flag indicating the 1st sample file in vcf format to be compared*  
-b | *flag indicating the 2nd sample file in vcf format to be compared*  
##### File Output 
The bedtools intersect tool outputs a vcf file of unique genome features based on the base file provided in the -a parameter.

## Python Script  

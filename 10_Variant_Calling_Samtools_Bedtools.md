# Variant Calling using *samtools* and *bedtools*  

## Commands and Parameters  
### samtools   
#### ***Library Index***
```
samtools faidx [reference file in fasta format]
``` 
#### ***mpileup***  
```
samtools mpileup -vf [reference file in fasta format] [sample file in sorted bam format] | bcftools call -vmO z -o [output file in vcf gzipped format]
```  
Tool | Command | Parameter | Description  
----|---------|-----------|------------  
samtools | mpileup | -v | *variant calling format (vcf).*  
| | | -f | *flag indicating reference file in fasta format.*  
bcftools | call | -v |
| | | -m |   
| | | -O |   
| | | z |   
| | |-o | *flag indicating output file path and file format.*  

#### File Ouput  
##### vcf.gz
### bedtools  
#### ***stats***  
```
bcftools stats --threads 11 -F [reference file in fasta format] -s [sample file in sam format] > [output file and path]
```  
Parameter | Description
----------|------------  
-threads | *number of threads available for use.*  
-F | *flag indicating reference file in fasta format.*  
-s | *flag indicating sample file in sam format.*  

##### ***plot-vcfstats***  
```
bcftools/plot-vcfstats -p [output folder name] [sample file from *bcftools stats* output file]
```  
Parameter | Description
----------|------------
-p |  
#### File Ouput 

# Python Script


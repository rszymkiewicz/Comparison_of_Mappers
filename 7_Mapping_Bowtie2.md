# Mapping using Bowtie2
Bowtie2 mapping and analysis will be performed on Sample P750. A python script will be used to perform bowtie2 mapping on all trimmed samples of the data set. The python script can be found at the bottom of the page.  

## Commands and Parameters  
### *Library Preparation*
```  
bowtie2-build [reference file in fasta format] [prefix of library] 
```  
### *Mapping*
```
bowtie2 -t -p 10 -x [prefix of library] -1 [Forward Paired fastq file] -2 [Reverse Paired fastq file] -S [Sam file output] --no-unal --no-discordant --no-mixed
```  
Parameter | Description  
----------|------------
-t | *flag indicating the computational power available.*  
-p | *flag indicating a paired-end factor.*  
-x | *flag indicating the previously structed library index.*    
-1 | *flag indicating the forward strand as an input file in fastq format.*  
-2 | *flag indicating the reverse strand as an input file in fastq format.*  
-S | *flag indicating the output file in sam format.*  
--no-unal |  
--no-discordant |  
--no-mixed |  

## Bowtie2 Analysis - Sample P750  

***Bowtie2 Sam Output***



# Python Script

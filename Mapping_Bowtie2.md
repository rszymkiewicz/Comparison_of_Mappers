# Mapping using Bowtie2
## Commands and Parameters  
### ***Library Preparation***
```  
bowtie2-build [reference file in fasta format] [prefix of library] 
```  
### ***Mapping***
```
bowtie2 -t -p 10 -x [prefix of library] -1 [Forward Paired fastq file] -2 [Reverse Paired fastq file] -S [Sam file output] --no-unal --no-discordant --no-mixed
```  
Parameter | Description |
-t | *number of threads* |
-p | *paired end* |
-x | *library* |
-1 | *forward strand in fastq format* |
-2 | *reverse strand in fastq format* |
-S | *file output in sam format* |
--no-unal | ** |
--no-discordant | ** |
--no-mixed | ** |

## Python Script

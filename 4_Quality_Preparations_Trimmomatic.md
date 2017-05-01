 #### Return to initial FastQC assessment of raw reads by clicking [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/3_Quality_Preparations_Inital_FastQC.md).  

# Quality Preparations - Trimming of Adapters and Low Quality Sequences

After assessing the inital quality of our sequencing samples we must take the appropriate steps to improve the quality of each sample to aid in sequence analysis. We will now quality trim our samples using a bioinformatic tool called ***Trimmomatic***.

## Trimmomatic 
**Trimmomatic**:  
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

**Trimmomatic** is a command line read trimming tool developed by the Usadel lab for both single-end and paired-end sequence data. Our sample set contains paired-end sequences indicated by R1 (forward strand) and R2 (reverse stand). In our case, we will use the paired-end features for trimming our sample set. A python script will be used to run trimmomatic on all samples of our data set and can be found at the bottom of the page.


Trimmomatic has two methods to to trim sequence reads:  
1.  ***Simple Trimming*** 
- The adaptor sequence is compared to the read and if it matches the minimum match threshold it is removed from the read.
2.  ***Palindrome Trimming***
- Designed for 'reading-through' a short fragment into the adaptor sequence on the other end. 
- The forward read is clipped of the adaptor sequence however the reverse read is dropped as it provides no further information. 

Trimmomatic requires two input fastq files (forward and reverse per sample) when performed paired-end trimming. There are a total of 4 output files produced by the trimmomatic paired-end feature. There are a total of two paired-end output files which contain the surviving paired reads output for each strand independently as well as, two unpaired output files which contain the reads which were unable to match to the subsequent strand but survived the trimming process.  

The command contains both adapter timming components and quality trimming components. The adaptor trimming component is defined by the ILLUMINACLIP parameter whereas the quality trimming components consist of the LEADING, TRAILING, SLIDINGWINDOW, and MINLEN parameters as described below.  

### Trimmomatic Code and Parameters
```javascript
java -jar trimmomatic.jar PE -threads 12 -trimlog [.txt] [R1.fastq][R2.fastq] [forward_paired.fastq.gz][forward_unpaired.fastq.gz] [reverse_paired.fastq.gz][reverse_unpaired.fastq.gz] ILLUMINACLIP:/home/rachelle/bin/trimmomatic-master/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```    


Parameter | Description  
----------|------------  
PE | *flag indicating Paired-end analysis*
-threads | *total number of computational power accessible for the software program to use.*  
-trimlog | *flag indicating the output text file containing a detailed log of all components trimmed.*  
-phred33 or -phred64 | *flag indicating the conversion of quality output to phred33 or phred64. Now performed automatically.*
**Adaptor Trimming**
ILLUMINACLIP | *flag indicating the fasta file containing the adapters to be removed.*  
ILLUMINACLIP:**2**:30:10| *The **first number** is the seedMismatches parameter which indicates the maximum mismatch count possible to be included upond trimming.*  
ILLUMINACLIP:2:**30**:10| *The **second number** is the palindromeClipThreshold parameter which indicates the threshold match required paired-end feature regarding two adapter ligated reads for the palindrome read alignment.*  
ILLUMINACLIP:2:30:**10**| *The **third number** is the simpleClipThreshold which indicates the minimum match theshold required to between the adapter sequence and the read being assessed.*  
**Quality Trimming**
LEADING | *removes leading low quality or N bases below the quality threshold indicated by the number (3).*  
TRAILING | *removes trailing low quality or N bases below the quality threshold indicated by the number (3).*  
SLIDINGWINDOW:**4**:15 |  *The size (number of bases) of the sliding window is indicated by the first number (4).*    
SLIDINGWINDOW:4:**15** |  *removes bases where the average quality of that base is less than the second number indicated (15).*
MINLEN | *removes reads which fall below the minimum number of bases per read (read length) threshold indicated by the number (3).*  

### Trimmomatic Analysis - Sample P750
The trimmomatic code above was used to perform adapter and quality trimming of Sample P750. The code states that we chose to remove Nextera adapter sequences, bases which had quality less than 3, those bases where the average quality of 4-bases wide was less than 15, as well as any reads less than 36 bases in length.  

#### *Trimmomatic Log Output*
The trimmomatic log output text file contains the following components: read name, suriving sequence length, location of first surviving base, location of last surviving base, and the amount trimmed from the end. All components are separated by a colon. To view the log output file created use the following code ```more [log_output_file name]``` to view the entire log output file or ```head [log_output_file_name]``` where the first 10 lines of the file are available as screen output.   

The first 10 lines of the log output file for Sample P750 is outlined below. We obtained this output by using the command ```head 
LOG_P750_raw_phred33```
```
==> LOG_P750_raw_phred33_2017April26 <==
HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 1:N:0:GCTACGCTCTAAGCCT 56 0 56 95
HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 2:N:0:GCTACGCTCTAAGCCT 151 0 151 0
HWI-1KL153:71:H9EA6ADXX:2:1101:3008:1954 1:N:0:GCTACGCTCTAAGCCT 123 0 123 28
HWI-1KL153:71:H9EA6ADXX:2:1101:3008:1954 2:N:0:GCTACGCTCTAAGCCT 0 0 0 0
HWI-1KL153:71:H9EA6ADXX:2:1101:5208:1883 1:N:0:GCTACGCTCTAAGCCT 150 1 151 0
HWI-1KL153:71:H9EA6ADXX:2:1101:5208:1883 2:N:0:GCTACGCTCTAAGCCT 151 0 151 0
HWI-1KL153:71:H9EA6ADXX:2:1101:5020:1893 1:N:0:GCTACGCTCTAAGCCT 150 1 151 0
HWI-1KL153:71:H9EA6ADXX:2:1101:5020:1893 2:N:0:GCTACGCTCTAAGCCT 139 0 139 12
HWI-1KL153:71:H9EA6ADXX:2:1101:5523:1877 1:N:0:GCTACGCTCTAAGCCT 150 1 151 0
HWI-1KL153:71:H9EA6ADXX:2:1101:5523:1877 2:N:0:GCTACGCTCTAAGCCT 0 0 0 0
```
The log output file states in line 1:  
**HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 1:N:0:GCTACGCTCTAAGCCT** 56 0 56 95 is the *sequence name*.  
HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 1:N:0:GCTACGCTCTAAGCCT **56** 0 56 95 is the *surviving sequence length*.  
HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 1:N:0:GCTACGCTCTAAGCCT 56 **0** 56 95 is the *location of the first surviving base*.  
HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 1:N:0:GCTACGCTCTAAGCCT 56 0 **56** 95 is the *location of the last suriviving base*.  
HWI-1KL153:71:H9EA6ADXX:2:1101:2421:1913 1:N:0:GCTACGCTCTAAGCCT 56 0 56 **95** is the *total amount trimmed from the end*.  

#### *Trimmomatic Screen Output*
Upon successful completion of trimmomatic on Sample P750, the following was output to the screen:  
```
Input Read Pairs: 1350885   
Both Surviving: 756928 (56.03%)   
Forward Only Surviving: 572214 (42.36%)   
Reverse Only Surviving: 8571 (0.63%)   
Dropped: 13172 (0.98%)  
TrimmomaticPE: Completed successfully  
```  
An explaination of the summary output for Sample P750 can be found within the chart below.  

Category | Explaination of Result  
---------|------------------------  
Input read pairs | *A total of **1,350,885 read pairs** were found between the forward and reverse strands used as input.*  
Both surviving | *A total of **756,928 reads pairs** survived the trimming process based on the outlined thresholds stated within the command above.*  
Forward Only Surviving | *A total of **572,214 reads from the forward strand** survived the trimming proceed set by the thresholds stated within the command above therefore; its paired read from the reverse strand failed one or more criteria of trimming.*  
Reverse Only Surviving | *A total of **8,571 reads from the reverse strand** survived the trimming proceed set by the thresholds stated within the command above therefore; its paired read from the forward strand failed one or more criteria of trimming.*  
Dropped | *A total of **13,172 paired reads** from both the forward and reverse strands failed to meet the trimming criteria therefore; they were both removed.*  

We can conclude that a total of 593,957 paired reads failed the thresholds set in the above command. In addition, less than 50% of the forward reads independently passed the thresholds and less than 1% of the reads from the reverse strand passed the threshold. Although this does seem like an extreme number lost in the trimming proccess we must keep in mind that the quality of our reverse reads were very low. Overall less than 1% of the paired reads were completely removed as they did not satisfy any of the independent or paired read conditions.  

# Python Script - Trimmomatic
## [Trimmomatic](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/Trimmomatic.py)

## Lets varify our quality and adapter trimming was successful and makes practical sense by using [FastQC](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/5_Quality_Preparations_FastQC_CheckAfterTrimming.md) again.  

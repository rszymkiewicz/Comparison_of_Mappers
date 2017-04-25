# Quality Preparations  

After assessing the inital quality of our sequencing samples we must take the appropriate steps to improve the quality of each sample to aid in sequence analysis. We will now quality trim our samples using a bioinformatic tool called ***Trimmomatic***.

## Trimmomatic 
**Trimmomatic**:  
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Trimmomatic is a command line read trimming tool developed by the Usadel lab for both single-end and paired-end sequence data. Our sample set contains paired-end sequences indicated by R1 (forward strand) and R2 (reverse stand). In our case, we will use the paired-end features for trimming our sample set. A python script will be used to run trimmomatic on all samples of our data set and can be found at the bottom of the page.

Trimmomatic requires two input fastq files (forward and reverse per sample). There are a total of 4 output files produced by the trimmomatic paired-end feature. There are a total of two paired-end output files which contain the surviving paired reads output for each strand independently as well as two unpaired output files which contain the reads which were unable to match to the subsequent strand but survived the trimming process.  

The command contains both adapter timming components and quality trimming components. The adaptor trimming component is defined by the ILLUMINACLIP parameter whereas the quality trimming components consist of the LEADING, TRAILING, SLIDINGWINDOW, and MINLEN parameters as described below.  

### Trimmomatic Code and Parameters
```javascript
java -jar trimmomatic.jar PE -threads 12 -trimlog [.txt] [R1.fastq][R2.fastq] [forward_paired.fastq.gz][forward_unpaired.fastq.gz] [reverse_paired.fastq.gz][reverse_unpaired.fastq.gz] ILLUMINACLIP:/home/rachelle/bin/trimmomatic-master/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:3
```    

Parameter | Description  
----------|------------  
PE | *flag indicating Paired-end analysis*
-threads | *total number of computational power accessible for the software program to use.*  
-trimlog | *flag indicating the output text file containing a detailed log of all components trimmed.*  
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

#### Trimmomatic Log Output
The trimmomatic log output text file contains the following components: read name, suriving sequence length, location of first surviving base, location of last surviving base, and the amount trimmed from the end. All components are separated by a colon.  

Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
Remove leading low quality or N bases (below quality 3) (LEADING:3) minimum quality required to keep a base (same for below)
Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
Drop reads below the 36 bases long (MINLEN:36)

# Python Script
## Trimmomatic 

# Lets varify our quality and adapter trimming was successful and makes practical sense by using FastQC again.  
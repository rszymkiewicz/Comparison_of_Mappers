# Mapping using Bowtie2  
***End-to-end Alignment and Paired-end***  
Bowtie2 performs alignment where all characters of a read must align to a segment of reference sequence. This process is known as end-to-end read alignment and is the default setting for Bowtie2. As our sample set is paired-end, strands of forward and reverse orientation, we will utilize the paired-end parameters indicated by -1 and -2 to represent the forward and reverse strand respectively. Bowtie2 has a randomized process for assigning alignment therefore when there are multiple possible alignments of good quality a pseudo-random interger determines which alignment will be reported.  

***Alignment Score***  
The alignment of each read to the reference sequence results in a score. The alignment score identifies the similarity between the sample sequence and the reference sequence. The end-to-end alignment score is determined by a series of points and penalties including that of a mismatched base or gap. The highest alignment score for the end-to-end alignment is 0 meaning the sample and the reference are the same with no differences between them. The alignment score must meet or exceed the minimum score threshold in order to be classified as valid alignment in Bowtie2. The threshold is calculated as a function of the read length ```(-0.6 + -0.6 x *read length*)```. The paired-end alignment score is the sum of all alignment scores of the individual strands of the pair, the forward and reverse strands.  

When an alignment is declared by Bowtie2, it is stating that the alignment is valid. When multiple alignments are declared by Bowtie2 they are said to be valid and distinct from each other. *Distinct* is defined as an individual read mapping (aligning) to different regions within the reference sequence. Bowtie2 looks for both distinct and valid alignments for all reads within its default settings. Bowtie2 will first search for valid alignments however; it will continue to look for alignments until it exceeds a search threshold or it cannot find another alignment which exceeds the current alignment score, from this MAPQ is calculated.  

It is important to note that the objective of Bowtie2 was computational efficiency not to output all possible data therefore it was not designed to report all alignments possible. If this is a desired outcome it is possible to perform this in Bowtie2 using the *-a* parameter which will output all possible alignments with the largest alignment score first.  

***Mapping Quality***  
Mapping quality is the aligner tools method to report the its confidence in assigning the read to its point of origin. As this is not always an easy task, as in situations of repeatitive regions where there is no preference to align in one segment versus another, the aligner is able to provide an estimation corresponding to its belief in its true point of origin. The mapping quality is abreviated to **MAPQ** and is output within the output *sam file*. Accurate mapping is essential to variant calling as it assists in identifying the uniqueness of the alignment between the sample and the reference sequence at specific regions across the sample. Low mapping quality regions can be ignored in variant calling as there maybe a potential that the read originated from another location rather than the one in question.  

***Concordant versus Discordant pairs***   
*Concordant alignment of pairs* is when the two strands of one pair align within the expected range of distance between the two strands. There are also some unique cases where Bowtie2 considers concordant alignment of pairs to be true if the two strands are overlapping or if they contain components of one another.  

If each strand as a unique alignment and the pairs do not align they represent a *discordant alignment of pairs.* When the pairs do not align this means that their orientations and/or expected distance range is not correct. These *discordant alignment of pairs* are important in structural variation as they have the potential to represent *inversions, translocations, insertions, or deletions*. By default Bowtie2 searches for both concordant and discordant alignments.  

Bowtie2 mapping and analysis will be performed on Sample P750. A python script will be used to perform bowtie2 mapping on all trimmed samples of the data set. The python script can be found at the bottom of the page.

## Commands and Parameters  
### *Library Preparation*
```  
bowtie2-build [reference file in fasta format] [prefix of library] 
```  

```bowtie2-build``` is a wrapper script which has the ability to index reference genomes of various sizes. Small indexes are build with *.bt2 file extension* whereas large indexes are build with *.bt21 file extension*.

### *Mapping*
```
bowtie2 -t -p 10 -x [prefix of library] -1 [Forward Paired fastq file] -2 [Reverse Paired fastq file] -S [Sam file output] --phred33 --no-unal --no-discordant --no-mixed
 ```  
Parameter | Description  
----------|------------
-t | *flag representing the timestamp*
-p | *flag indicating the computational power available. Allows for parallel processing and high-throughput.*   
-x | *flag indicating the previously structed library index.*    
-1 | *flag indicating the forward strand as an input file in fastq format.*  
-2 | *flag indicating the reverse strand as an input file in fastq format.*  
-S | *flag indicating the output file in sam format.*  
--phred33 | *flag indicating sample file inputs are of the latest Illumina pipelines at phred quality 33.*
--no-unal | *flag utilized for sam file output option; not to include reads which did not align to reference sequence.*  
--no-mixed | *flag indicating that Bowtie2 will only search for paired-end alignment and not used unpaired alignments from individual strands.*  

***Default Bowtie2 settings have *sensitive* preset parameters which allows for more accurate and sensitive results.***  
```--sensitive [-D 15 -R 2 -L 22 -i S,1,1.15]```  

Default Parameter | Description  
------------------|------------
-D | *flag indicating the number of seed extension attempts to identify the best alignment of a read to the reference sequence ('total number of failures') prior to obtaining a new read to align. A failure is defined as being unable to exceed the current best alignment score.*  
-R | *flag indicating re-seeding reads with repetitive seeds. Reads have repetitive seeds when the total number of seed hits divided by the number of seeds which aligned at least once is greater than 300.*  
-L | *flag indicating the length of seed substring for multiseed alignment.*  
-i | *flag indicating the interval between seed substrings for multiseed alignment.*  

## Bowtie2 Analysis - Sample P750
***Bowtie2 Alignment Summary***
```Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 00:01:03
756928 reads; of these:
  756928 (100.00%) were paired; of these:
    87782 (11.60%) aligned concordantly 0 times
    651108 (86.02%) aligned concordantly exactly 1 time
    18038 (2.38%) aligned concordantly >1 times
    ----
    87782 pairs aligned concordantly 0 times; of these:
      38250 (43.57%) aligned discordantly 1 time
93.46% overall alignment rate
Time searching: 00:01:04
Overall time: 00:01:04
```  

Mapping of the 756,928 reads of the Sample P750 to the reference sequence using Bowtie2 (see above command) produced an **overall alignment of 93.46%**. Bowtie2 alignment is broken down into three categories: concordant, discordant, and those reads which did not align at all or aligned indivdiually. In the P750 Sample, there are a total of 669,146 reads (651,108 + 18,038) which align to the reference concordantly (88.40%). Of the reads which did not align in a concordant manner (87782 reads), 38250 of the 87782 reads aligned discordantly (43.5%). The remaining reads are the reads which did not align as paired-ends therefore; they aligned as single reads or may not have aligned at all to the reference. A total of 49,532 reads fall into this category: total - concordant - discordant (756,928 - 669,146 - 38250). To determine the overall alignment because this is paired-end alignment, we must multiple our reads by 2 to represent each strand. The overall alignment is calculated by twice the number of concordant aligned reads plus twice the number of discordant number of aligned reads divided twice the total number of reads from the sample (((2 x 669146) + (2 x 38250))/1513856).

***Bowtie2 Sam Output***  
The *sam file* for paired-end alignment has information for each strand of the pair on an individiual line. The first line represensts the alignment for the forward strand whereas the second line represents the alignment of the reverse strand. The *sam file* has multiple fields however; the RNEXT and PNEXT fields contain information pertaining to the paired-end alignment specifically the reference name and position to which the other strand aligned to the reference sequence. 

Field # | Description  
---------|------------  
2 | *describes the paired-end nature of read and alignment*

# Python Script - Bowtie2 Mapping

## Now lets try using [BWA](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/8_Mapping_BWA.md) to map our samples to our reference sequence.

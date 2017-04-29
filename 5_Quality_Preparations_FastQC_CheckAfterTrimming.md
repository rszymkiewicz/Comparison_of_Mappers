# Quality Preparations - FastQC Quality Check After Trimming
Upon the completion of adapter and quality trimming of raw sequence reads it is important to verify trimming was not only successful but accurate in its performance. In some instances, trimming may not have resulted in the desired resultant. By utilizing FastQC software program one is able to confirm the success or further quality preparations required prior to continuing with sequence analysis. A python script to perform post trimming sequence quality check for the entire data set can be found at the bottom of the page. We will continue to follow our P750 sample and assess the success of our previous trimming step using trimmomatic below.  

**FastQC:**  
Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc  

## FastQC Quality Check After Trimming Command
```
fastqc [forward_trimmed.fastq] -o [output_folder_path]  
fastqc [reverse_trimmed.fastq] -o [output_folder_path]  
```  

Parameter | Description  
----------|------------  
-o | *flag indicating output folder path location*  

## FastQC Analysis - Sample P750  
As previously explored in our initial FastQC assessment of our raw samples from our sample set, we focused on 3 specific areas of analysis; Basic Statistics, Per base sequence quality, and Adapter Content. Post-trimming we will look specifically at these 3 areas again for our Sample P750 however; we will also take into account any other areas which may have been flagged as unusual (red) and draw conclusions surrounding our trimming step. A total of 3 of the 12 areas of analysis were marked as red (unusual), *per tile sequence quality*, *per base sequence content*, and *Kmer content*. We will also assess these 3 areas of analysis and determine if we can proceed forward with our sample set.  

### Basic Statistics  
pIt is important to verify no extreme results were obtained after trimming of our sample. 
![Basic_Statistics_Post_Trimming_R1](https://cloud.githubusercontent.com/assets/25803304/25551963/e9075da0-2c5b-11e7-9f3b-1795c0309557.png)  
![Basic_Statistics_Post_Trimming_R2](https://cloud.githubusercontent.com/assets/25803304/25551958/e8f63192-2c5b-11e7-92d7-461c728e9ce2.png)  

***Analysis***  
It is evident there will be a decrease in the number of reads in both our forward and reverse strand. For our P750 sample, the forward  and reverse strand lost a total 593,957 reads. Our sequence length range for each raw read has become a range of 36 to 151bp for both the forward and reverse strands. Post-trimming the %GC content for both forward and reverse strands was identified to be 65.  

### Per Base Sequence Quality  
![Per_Base_Sequence_Quality_Post_Trimming_R1](https://cloud.githubusercontent.com/assets/25803304/25551962/e9065734-2c5b-11e7-8192-2e31a324a4e2.png)  
![Per_Base_Sequence_Quality_Post_Trimming_R2](https://cloud.githubusercontent.com/assets/25803304/25551955/e8f2612a-2c5b-11e7-8d28-303878201e11.png)  

***Analysis***  
All reads of both the forward and reverse strands appear within the green region of the graph therefore; deemed to be high quality reads.  

### Adapter Content  
Here we verify that the adapter sequences seen in the inital FastQC output prior to performing trimming are now gone. If this is not the case then we must continue to assess and remove unwanted library sequences.  
![Adapter_Content_Post_Trimming_R1](https://cloud.githubusercontent.com/assets/25803304/25551960/e90106b2-2c5b-11e7-91f2-fab911e8e456.png)  
![Adapter_Content_Post_Trimming_R2](https://cloud.githubusercontent.com/assets/25803304/25557770/316fd740-2ce6-11e7-9b95-a764119f154d.png)  

***Analysis***  
In the case of Sample P750, all Nextera XT Illumina Library Kit sequences have been successfully removed.  

### Per title Sequence Quality
![Per Title Sequence Quality_Post_Trimming_R1](https://cloud.githubusercontent.com/assets/25803304/25551965/e90e2356-2c5b-11e7-98bd-f4dc7171f4c1.png)  
![Per Title Sequence Quality_Post_Trimming_R2](https://cloud.githubusercontent.com/assets/25803304/25551957/e8f5dfb2-2c5b-11e7-85d8-5f8ec01fae2a.png)  

***Analysis***  
Trimmomatic manual states that the **per title sequencing quality** calls a failure when the phred score is 5 less than the mean for that base across all titles. Reasons for seeing warnings or errors on this plot could be transient problems such as bubbles going through the flowcell, or they could be more permanent problems such as smudges on the flowcell or debris inside the flowcell lane. The forward strand has a few areas of moderate/high concern regarding a low phred quality score however; overall there is no concern. The reverse strand on the other hand has a large number of titles which low phred quality score. This may be attributed by the low quality reads which are common in reverse strands.  

### Per Base Sequence content  
Demonstrates the proportion of each of the four DNA bases been called for each base position in a sample file.  
![Per Base Sequence Content_Post_Trimming_R1](https://cloud.githubusercontent.com/assets/25803304/25551961/e9063ee8-2c5b-11e7-8fd6-1ef255605c42.png)  
![Per Base Sequence Content_Post_Trimming_R2](https://cloud.githubusercontent.com/assets/25803304/25551956/e8f524f0-2c5b-11e7-8b5d-be00ed7c30da.png)   

***Analysis***  
As FastQC assumes a random library there should be minimal if any bias in nucleotide calling. FastQC manual states that **per base sequence content** calls a failure when a difference between A and T, or G and C is greater than 20% in any position. In some instances, certain library types will always produce biased sequence composition, normally at the start of the read and that in most cases it does not affect downstream analysis. Here we do see a biased sequence composition at the start of the read for both the forward and reverse strand. It appears that the first 15bp of the sample reads are biased. The FastQC manual also states there is minimal impact to downstream analysis but as the bias is quite high we will need to keep this in mind during analysis.  

### Kmer Content  
![Kmer Content Post Trimming R1](https://cloud.githubusercontent.com/assets/25803304/25551954/e8e72a26-2c5b-11e7-878f-7974f91617b9.png)  
![Kmer Content Post Trimming R2](https://cloud.githubusercontent.com/assets/25803304/25551959/e8f78470-2c5b-11e7-8042-4f5da95b5c5b.png)  

***Analysis***  
Trimmomatic manual states that the **kmer content** will fail if the k-mer is imbalanced with a binomial p-value <10^5. It also states that libraries which derive from random priming will nearly always show Kmer bias at the start of the library due to an incomplete sampling of the possible random primers. The information from the kmer content will be kept in mind during downstream analysis.  

## Summary Analysis - Sample P750:  
The Nextera XT DNA Library adapter sequences were successfully removed in trimming and low quality reads were also removed in trimming. Although 3 of the 12 areas of analysis were marked as red (unusual), *per tile sequence quality*, *per base sequence content*, and *Kmer content* the warning threshold is deemed to be sensitive for these 3 categories as described within the FastQC manual. We will keep these warnings in mind throughout the remaining downstream analysis however; we will proceed forward in confidence with our trimmed P750 sample.  

# Python Script - FastQC Quality Check Post Trimming

## Now that we have confidence in our sample sequences we can continue moving forward with our analysis. Lets look further into performing [*mapping*](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/6_Mapping.md) of our samples to our reference sequence.   

# Quality Preparation of Whole Genome Sequences

Checking for good quality sequences is the first component of sequence analysis.  
Quality checking can be performed through the evaluation of raw sequences. To do this, programs can be used to identify potential areas of uncertainty specifically surrounding both the sequencer (equipment error) and the library material (human error) which may have contributed to the final error rate. It is essential to investigate the quality of sequences, as the biological outcomes will be affected by each step of analysis, beginning with the raw sample. 

In order to begin sequence analysis with the highest sample quality possible, one must first assess the sequences followed by trimming of inappropriate components such as library material and low quality reads. The FastQC software tool is a great starting point to identify potential areas of uncertainity within the sequence sample.  

## FastQC
**FastQC:**  
Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

FastQC a java application developed by Babraham Bioinformatics which identifies and categorizes potential problems with raw sequence data from high throughput sequencing pipelines such as Illumina. The FastQC application provides an easy to understand html summary report covering 12 areas of interest. Each area of analysis is evaluated and measured based on the assumption of a random and diverse sample. A pass (normal), warning (suspicious), and failure (unusual result) of analysis indicated by a green, yellow, or red colour respectively. It is important to note that all areas of analysis should be taken in the context of the specific sample and its conditions. We will be using this program to evaluate the raw reads in fastq format of our sample set focusing on sample P750. A python script used to perform FastQC inital quality check on all samples within the sample set can be seen at the end of this page. 

### FastQC Quality Check Command 
```
fastqc [Forward_Strand.fastq] -o [output_folder_path]
fastqc [Reverse_Strand.fastq] -o [output_folder_path]
```  
Parameter | Description  
----------|-------------  
-o | *flag indicating the output folder path which will contain all output files including the html file*  

### FastQC Analysis - Sample P750
Looking at our specific sample (Sample P750) we can focus on specific areas to evaluate in order to assess aspects which require quality trimming. Here we will look at 3 specific areas of analysis; *Basic Statistics*, *Per base sequence quality*, and *Adapter Content*.  

#### Basic Statistics
*Basic Statistics* is a summary in chart format of the details and important information regarding the analyzed sample. Specific attention should be brought to the *encoding* aspect of the chart as it incidates the ASCII encoding of quality values for the analyzed sample.  

![Basic_Statistics_Forward_Strand_R1](/media/sdc/P750_Github/P750_R2_fastqc_before.png)  
![Basic_Statistics_Reverse_Strand_R2](/media/sdc/P750_Github/P750_R2_fastqc_before.png)  

***Analysis***  
From the *Basic Statistics* chart we can conclude that there are a total of 1,350,885 raw reads each within the forward and reverse strands in addition to the sequence length of each raw read was 151bp in length. The average %GC content for the forward strand was slightly higher than the reverse strand (65% versus 64%) which is suitable range for *Pseudomonas aeruginosa* strains. In addition, the encoding parameter tells us that the sequencing is phred33 encoding for both forward and reverse strands.   

#### Per base sequence quality
*Per base sequence quality* is a boxplot statistical summary of quality values across all bases at each position in the input fastq sample file. The y-axis represents the quality scores and the x-axis represents the base position. A key aspect of this graph is the 3-colour divisions which help guide the interpretation of results. The green background indicates good quality calls where the red background indicates poor quality calls. It is important to note that the quality degrades as near the end of the sequencing run therefore in most cases it is unlikely to see high quality reads at the end of the read.  

![Per_base_sequence_quality_Forward_Strand_R1](/media/sdc/P750_Github/P750_R1_perqualitybase_fastqc_before.png)  
![Per_base_sequence_quality_Reverse_Strand_R2](/media/sdc/P750_Github/P750_R2_perqualitybase_fastqc_before.png)  

***Analysis***  
*Forward Strand (R1)*  
The forward strand *per base sequence quality* was successful as a green check mark appears for this category. It can be seen within the figure that almost all reads are found to have high quality scores across all base positions as indicated by the boxplots laying within the green region.  

*Reverse Strand (R2)* 
The reverse strand *per base sequence quality* however; was not successful as a red 'x' mark appears for this category. It can be seen within the figure that almost all reads are found to have medium to very low quality scores across all base positions as indicated by a majority of the boxplots laying within the orange and red regions.  

#### Adapter Content
*Adaptor Content* is a linear graph presenting the adapter sequences found present within the input sample. The FastQC tool 
has the ability to detect 4 types of adapter sequences including the Nextera Transposase Sequence. In our sample set, the sequencing facility utilizes a Nextera Library Kit therefore it is possible to see Nextera adapter sequences post sequencing. Adapter libraries should be removed in quality analysis.  

![Adapter_Content_Forward_Strand_R1](/media/sdc/P750_Github/P750_R1_adaptercontent_fastqc_before.png)
![Adapter_Content_Reverse_Strand_R2](/media/sdc/P750_Github/P750_R2_adaptercontent_fastqc_before.png)

***Analysis***  
Both the forward and reverse strands have presence of adapter sequences. As the Nextera XT was used for sequencing this sample set it is natural to see Nextera sequencing appear in the adapter content section of the FastQC analysis. This adapter sequence must be removed prior to completing downstream analysis.   

# Python Scripts
## FastQC Quality Check - Python Script

# Now that we have assessed our raw sample reads lets perform the appropriate trimming of our samples using [Trimmomatic](https://github.com/rszymkiewicz/Comparison_of_Mappers/edit/master/4_Quality_Preparations_Trimmomatic.md).


# Quality Preparation of Whole Genome Sequences

Checking for good quality sequences is the first component of sequence analysis. Quality checking can be performed through the evaluation of raw sequences. To do this, programs can identify potential areas of uncertainty specifically surrounding both the sequencer (equipment error) and the library material (human error) which may have contributed to the final error rate. It is essential to investigate the quality of the sequences prior to analysis as the biological outcomes will be affected by each step of analysis beginning with the raw sample. 

In order to begin sequence analysis with the highest sample quality possible, one must first assess the sequences and performing trimming of inappropriate components such as library material and low quality reads. The FastQC software tool is a great starting tool to identify potential areas of uncertainity within the sequence sample.  

## FastQC
Is a java application which identifies categories of potential problems with raw sequence data from high throughput sequencing pipelines such as Illumina. The FastQC application provides an easy to understand html summary report covering 12 areas of interest. The FastQC program has an easy to understand evaluation of results indicated by a green, yellow, or red colour per area of analysis and are measured based on the assumption of a random and diverse sample. Green indicates normal, yellow indicates suspicious results, and red indicates unusual results. All areas of analysis should be taken in the context of your sample and conditions. Specific areas of focus include the *Basic Statistics*, *Per base sequence quality*, and *Adapter Content*. We will be using this program to evaluate the raw reads in fastq format of our sample set.    

### Basic Statistics
*Basic Statistics* is a summary in chart format of the details and important information regarding the analyzed sample. 
![Basic_Statistics](Full path of image)  
Specific attention should be brought to the *encoding* aspect of the chart as it incidates the ASCII encoding of quality values for the analyzed sample.  

### Per base sequence quality
*Per base sequence quality* is a boxplot statistical summary of quality values across all bases at each position in the input fastq sample file. The y-axis represents the quality scores and the x-axis represents the base position. A key aspect of this graph is the 3-colour divisions which help guide the interpretation of results. The green background indicates good quality calls where the red background indicates poor quality calls. It is important to note that the quality degrades as near the end of the sequencing run therefore in most cases it is unlikely to see high quality reads at the end of the read.  
![Per_base_sequence_quality](Full path of image)  

### Adapter Content
*Adaptor Content* is a linear graph presenting the adapter sequences found present within the input sample. The FastQC tool 
has the ability to detect 4 types of adapter sequences including the Nextera Transposase Sequence. In our sample set, the sequencing facility utilizes a Nextera Library Kit therefore it is possible to see Nextera adapter sequences post sequencing. Adapter libraries should be removed in quality analysis.  
![Adapter_Content](Full path of image)

### Quality Check Command 
```
fastqc [Forward_Strand.fastq] -o [output_folder_path]
fastqc [Reverse_Strand.fastq] -o [output_folder_path]
```  
### Quality Check Results:

### Quality Check Post-Trimming Results:

## Trimmomatic

### Trimmomatic Command

# Python Scripts
## Quality Check
## Trimming using Trimmomatic

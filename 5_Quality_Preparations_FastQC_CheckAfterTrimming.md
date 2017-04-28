# Quality Preparations - FastQC Quality Check After Trimming

Upon the completion of adapter and quality trimming of raw sequence reads it is important to verify trimming was not only successful but accurate in its performance. In some instances, trimming may not have resulted in the desired resultant. By utilizing FastQC software program one is able to confirm the success or further quality preparations required prior to continuing with sequence analysis. A python script to perform post trimming sequence quality check for the entire data set can be found at the bottom of the page. We will continue to follow our P750 sample and assess the success of our previous trimming step using trimmomatic below.  

## FastQC Quality Check After Trimming Command
```
fastqc [forward_trimmed.fastq] -o [output_folder_path]  
fastqc [reverse_trimmed.fastq] -o [output_folder_path]  
```  

Parameter | Description  
----------|------------  
-o | *flag indicating output folder path location*  

## FastQC Analysis - Sample P750  
### Basic Statistics  
***Analysis***  
![Basic_Statistics_POST] (Full path to image)  

### Per Base Sequence Quality  
***Analysis***  
![Per_Base_Sequence_Quality_POST] (Full path to image)  

### Adapter Content  
***Analysis***
Here we verify that the adapter sequences seen in the inital FastQC output prior to performing trimming are now gone. If this is not the case then we must continue to assess and remove unwanted library sequences.  
![Adapter_Content_POST] (Full path to image)  

In the case of Sample P750, all Nextera XT Illumina Library Kit sequences have been successfully removed.  

# Python Script - FastQC Quality Check Post Trimming

## Now that we have confidence in our sample sequences we can continue moving forward with our analysis. Lets look further into performing [*mapping*](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/6_Mapping.md) of our samples to our reference sequence.   

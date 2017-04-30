# Mapping Output - Creating a sorted bam file

The ***sam*** file output which was created from each of the two mapping software tools [Bowtie2](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/7_Mapping_Bowtie2.md) and [BWA](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/8_Mapping_BWA.md) is a human readable output file however; it is quite large in size. Therefore in order to reduce the amount of memory being used by the computer, we can convert the *sam file* to its binary file counterpart known as a ***bam*** file. Once a *bam* file is created it can also be sorted by alignment position. A ***sorted bam*** file is an important aspect which will be used for downstream analysis such as variant calling.  

To do this we will use **```samtools view```** command in order to convert our *sam* output files from our mapping steps into a *bam* file. We will then sort and index our sorted *bam* file for quicker random access to regions of particular read alignments. 

## Command and Parameters
***samtools view***  
samtools view converts the human readable *sam* file into a computer readable (more efficient) binary file known as a *bam* file.  

```
samtools view -S -h -F 4 -b [sam_file] > [bam_file]
```

Parameter | Description  
----------|------------  
-S | *flag indicating to ignore any incompatibilities between samtools versions. This is a required flag dependent on the version of samtools. Later versions of samtools identify incompatibilities automatically and make the necessary modifications.*  
-h | *flag indicating to include sam header in coverstion.*  
-F | *flag indicating not to include alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].*  
-b | *flag which indicates the conversion of sam file format to bam file format.*  

### Note: upon creation of a **bam** file the **sam** file can be deleted as the binary file is more favourable computationally. If a *sam* file is needed at any point the *bam* file can be converted to a *sam* file.  
The following command can be used to remove the *sam* file from the directory:  
```
rm -rf [sam_file]
```  

***samtools sort***  
**samtools sort** sorts alignments based on their alignment position. The second parameter within the samtools sort command represents the prefix name of the output file to be created and the *.bam* extention is added to the end of the prefix defined by the user in the command.  

```
samtools sort [bam_file] + [sorted_bam_file_prefix]
```  
## We have completed alignment of our samples to our reference sequence. Lets continue and assess the possible SNPs and indels which maybe present within our sample set in comparison to our reference sequence. To do this lets learn about [*variant calling*](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/9_Variant_Calling.md).

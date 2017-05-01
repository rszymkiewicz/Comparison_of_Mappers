#### Return to *Data Comparison using vcftools and bedtools* by clicking [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/13_Comparison_vcftools_bedtools.md).  

# Comparison of Mappers Tutorial Conclusion

We have completed the variant calling and comparison of outputs for our two mapping tools, Bowtie2 and BWA.  
Lets review our findings.

## Mapping Tool Summary
Bowtie2 | BWA-MEM  
--------|----  
**93.46%** overall alignment | **99.08%** overall alignment  

We would expect to see similar overall alignment between the two mapping tools as they both are BWT algorithim based mapping tools. The approximately 5% difference in alignment could be a result of specific parameters and thresholds set by the user when preparing the commands. It is also important to note that the Bowtie2 command used in this tutorial called end-to-end alignment whereas traditionally BWA by default calls local alignment. Recall that end-to-end alignment requires all bases to align exactly whereas local does not require all bases to match rather a minimum threshold. Details of each tool can be found on their respective tutorial pages [Bowtie2](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/7_Mapping_Bowtie2.md) and [BWA-MEM](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/8_Mapping_BWA.md).  

## Variant Calling Comparison Summary
Comparison Tool | Mapping Tool | Results  
----------------|--------------|--------  
*vcf-compare* | **Bowtie2** | **28** *unique* calls  
*vcf-compare* | **BWA-MEM** | **14** *unique* calls  
*vcf-compare* | **Bowtie2 and BWA-MEM** | **77** *common* calls  
*bcftools subtract* | **Bowtie2** | **27** *unique* calls  
*bcftools subtract* | **BWA-MEM** | **14** *unique* calls  
*bcftools intersect* | **Bowtie2 and BWA-MEM** | **78** *common* calls  

Two of the comparison tools in the final stage of the tutorial produced very similar variant call results. Bowtie2 has a increased number of variant calls (approximately 14 calls) in comparison to BWA-MEM. Further analysis will need to take place to identify the biological significance and whether or not the variant calls made for each mapping tool are in fact real. 

## Biological Relevance
Sample P750 is a PES sample from the same patient as our reference sample. With this known, we expect to see similarity between Sample P750 sequence and the reference sequence. Sample P750 was taken 19 months after the reference sequence sample was taken. 

## Final remarks
As mentioned within the [project outline](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/1_Overview_of_Project.md) page of the tutorial, success of the project would be evaluated based on the following aspects:  
### **Mapping:**  
  A. Ability to obtain a sam output for each of the mapping tools.  
  B. Describe summary output of mapping tools.  

**Success:** *As this tutorial was able to perform mapping using two tools (```Bowtie2``` and ```BWA_MEM```) in order to produce a sam file output as well as, analyze the summary output provided by both mapping tools upon completion of alignment these objectives were met.*   

### **Variant Calling:**  
  A. Ability to obtain a vcf file for each of the mapping tools.  
  B. Describe features obtained in each of the vcf output files.  
  C. Identify the total number of SNPs and indels per mapping tool.  
  D. Identify similarities and differences in variant calling (SNP and indel calls) between the two mapping tools.  

**Success:** *As this tutorial was able to perform variant calling using ```samtools mpileup``` and ```bcftools call``` in order to prduce a vcf output file as well as identify features of a vcf file, the first two objectives of this section were met. Identifying the total number of SNPs and indels for each mapping tool was performed using ```bcftools stats``` and ```bcftools\plot-vcfstats``` and similarities and differences between variant calling was performed using ```vcf-compare``` and ```bedtools intersect/subtract```. Summary information was obtained for all comparisons therefore these project objectives were also met.*

In conclusion all desired project objectives were satisfied. Further analysis and comparisons pertaining to all samples of the PES sample set will be performed outside of this tutorial.  

## Acknowledgements
I would like to thank my professors for their guidance throughout the BIO720 course: Dr. Ian Dworkin, Dr. G. Brian Golding, and Dr. Ben J. Evans. 

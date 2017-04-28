# Project Objective:
Students are given the ability to design an independent project to analyize sequence data utilizing various bioinformatic software tools and suites learnt throughout the BIO720 course. Students are able to use any sequence data readily accessible for analysis.

I will create a high throughput bioinformatic pipeline that will incorporate both quality preparations and analysis components. The pipeline will be a combination of python scripts which will include a variety of bioinformatic tools and suites. **The objective of this pipeline is to compare two mapping/alignment tools (Bowtie2 and BWA) to identify similarities and differences in variant calling between the two mapping tools and to see which mapping tool provides the most biologically accurate SNPs and indels.**   

The pipeline consist of 4 components (Figure 1):
 - Quality Preparations
 - Mapping/Alignment
 - Variant Calling
 - Comparison of Results
 
 ![Figure1](C:/Users/dance/Pictures/Figure1_Github_BIO709CourseProject.png?raw=true)

Specifically, I will use this pipeline to assess 206 Prairie Epidemic Strains (PES) from cystic fibrosis patients to uncover SNPs and indels.

# Project Aims:
1. Gain experience with commonly used bioinformatic software tools used in analysis such as Bowtie2, BWA, samtools, bcftools, vcftools, and bedtools.
2. Compare results of mapping tools available for Whole Genome Sequence data. 
	A. Preparing samples for downstream analysis by performing quality steps such as trimming and removal of low quality reads from sequence data.
	B. Perform alignment of whole genome sequences to a reference sequence using Bowtie2 and BWA.
	C. Perform variant calling and compare vcf outputs of the two selected mapping tools. 

# Project Evaluation:  
Successful completion of this project will be determined based on the follow two major areas of the project:  
1. ***Mapping:***  
	A. Ability to obtain a vcf file for each of the mapping tools.  
	B. Describe features obtained in each of the vcf output files.   
2. ***Variant Calling:***  
	A. Identify the total number of SNPs and indels per mapping tool.  
	B. Identify similarities and differences in variant calling (SNP and indel calls) between the two mapping tools.  

## The project tutorial will begin with gaining some information about the selected Illumina sequencing [sample set](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/2_Sample_Background.md).

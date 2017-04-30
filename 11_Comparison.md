# Data Comparison

## Objective of this Project
***Goals of Comparison***  
Lets recap what we have done thus far. ***The objective of this pipeline as outlined in the [first tutorial page](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/1_Overview_of_Project.md) was to compare two mapping/alignment tools (Bowtie2 and BWA) by identifying similarities and differences in variant calling between the two mapping tools.***  

We have taken raw Illumina sequences and performed quality trimming, followed by mapping of our samples to our reference sequence. Using a sorted bam file we then performed variant calling for both mapping tools. We will now continue moving forward by comparing the results of our variant calling from our two mapping methods.  

## Tools Used for Comparison  
We will use a variety of bioinformatic tools to assess the similarities and differences in variant calling.  
1. vcftools - *vcf-compare*  
2. bcftools - *bcftools stats*  
3. bcftools - *bcftools/plot-vcfstats*  
4. bedtools - *bedtools intersect*  
5. bedtools - *bedtools subtract*  

## Now that we have combined all the identified the regions of alignment between our sample and reference sequence, lets dive deeper and compare the results of our two mapping tools (Bowtie2 and BWA) using [vcftools, bcftools, and bedtools](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/12_Comparison_vcftools_bedtools.md) to see if we obtained different results.  

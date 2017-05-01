#### Return to *Mapping Output File Conversion* by clicking [here](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/9_sam_to_bam_to_sorted.md).  

# Variant Calling
## What is *variant calling*?  
Variant calling of next-generation sequencing (NGS) data is a process to which raw sample sequence(s) are mapped to a reference sequence in order to identify similarities and differences that may occur between the raw sample sequence(s) and the reference sequence.  

A particular file known as the *variant calling format (VCF)* established by Danecek et al, in 2011 and is a file format commonly used to store variant calling information such as SNPS, indels, and structrual variants. VCFtools is a bioinformatic suite tool which manipulates this type of data file; allowing users to perform various analysis by merging or comparing *vcf* files. The manual for vcftools can be found [here](https://vcftools.github.io/index.html). 

## Citation
Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., … 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156–2158. http://doi.org/10.1093/bioinformatics/btr330

## Prior to manipulating a vcf file we first must create one. Lets use [samtools](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/11_Variant_Calling_Samtools_Bcftools.md) to help us do this.  

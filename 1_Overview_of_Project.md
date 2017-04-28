# Project Objective:
Students are given the ability to design an independent project to analyize sequence data untilizing various bioinformatic software tools and suites learnt throughout the course. Students are able to use sequence data readily accessible for analysis.

I will create a high throughput bioinformatic pipeline that will incorporate both quality preparations and analysis components. The pipeline will be a combination of python scripts including various bioinformatic tools and suites. The objective of this pipeline is to compare 2 mapping/alignment tools (bowtie2 and BWA) to identify which alignment will provide the most biologically accurate SNPs and indels. 
The pipeline consist of 4 components (Figure 1):
 - Quality Preparations
 - Mapping/Alignment
 - Variant Calling
 - Result Comparison (VCF Files)
 
 ![Figure1](C:/Users/dance/Pictures/Figure1_Github_BIO709CourseProject.png?raw=true)

Specifically, I will use this pipeline to assess 206 Prairie Epidemic Strains (PES) from cystic fibrosis patients to uncover SNPs and indels.

# Project Aims:
1. Gain experience with commonly used bioinformatic software tools such as bowtie2, BWA, samtools, and bedtools.
2. Compare results of mapping tools available for Whole Genome Sequence data. 
	A. Preparing samples for downstream analysis by performing quality steps such as trimming on sequence data.'
	B. Perform alignment of whole genome sequences to a reference sequence using bowtie2 and BWA.
	C. Perform variant calling and compare vcf outputs of 2 mapper tools. 
3. Project success will be evaluated based:
	A. Mapping effectiveness
		1. Identify which mapping tool had a higher alignment rate
	B. Variant Calling
		1. Total number of SNPs and indels identified by the mapper tool.
		2. Dissimilarities between SNP and indel calls between the two mapping tools.

# The project tutorial will begin with gaining some information about the selected Illumina sequencing [sample set](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/2_Sample_Background.md).

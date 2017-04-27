# General Information
## *What is PES?*
PES stands for *Prairie Epidemic Strain* and is a *Pseudomonas aeruginosa* bacterial strain found within Canadian cystic fibrosis adult patients at the Calgary Adult Cystic Fibrosis Clinic (CACFC). 

## *Cystic Fibrosis PES Samples*
The CACFC has maintained a biobank of cystic fibrosis sputum samples for over 50 years. The Surette laboratory has access to this unique biobank and performs research pertaining to various aspects of cystic fibrosis disease in adult patients. 
With the 206 PES samples available to the Surette laboratory we will be performing a longtiduinal study investigating the population genetics of PES over 30 years. 

# 206 PES Genome Sequences
## Illumina Sequencing
The 206 PES genome sequences are a combination of HiSeq and MiSeq Illumina sequences completed in 2013 and 2016 respectively. The Nextera XT DNA Library Prep Kit was used for both Illumina runs. An important aspect to note between these two Illumina sequence runs are the difference in read length; the Hiseq has on average 150bp raw read size while the MiSeq run consists of average read lengths of 251bp.  

## Metadata
There are a total of 206 PES genomes readily available to the Surette Laboratory for the population genetics longitudinal study. The 206 PES samples are divided into two groups, MDR [multi-drug resistant] (85 samples) and non-MDR (121 samples).  A total of 58 patients are represented across these two groups, 24+2* MDR patients and 56 patients for non-MDR samples. Each of the patients have a minimum of 2 samples (paired early and late isolates).  
***Note*** * *: 24 patients are shared across the two sample sets. 2 patients do not have non-MDR samples.* 

## Reference Sequence
P749 is the oldest sample of PES presently available in the Surette Laboratory. It date of origin is 1980. This P749 sample was sequenced by 2 different genome sequencers, Illumina Miseq and PacBio Sequencer. The two samples were prepared by Saad Syed for genomic sequencing. 

The P749 PES sample was selected to be sequenced by PacBio as 
PacBio sequencing is said to be the golden standard in genomic 
sequencing. The purpose of utilizing a reference sequence in 
analysis versus *de novo approach* is to compare unknown samples with the closest well defined (annotated) sample possible. In our case, the P749 PacBio sequenced sample will be utilized as our reference sequence throughout the analysis process of these 206 PES samples as it is at the highest standard of genomic sequencing. Details pertaining to the P749 PacBio Sequence can be seen in below.  
***P749 PacBio Sequence***
 - 64 216 632 bp
 - 165X Coverage
 - 3 contigs
 - 66.46% GC content

# Now that we know a bit more about our sample set and reference sequence lets begin by assessing the quality of our Illumina sample sequences using [FastQC](https://github.com/rszymkiewicz/Comparison_of_Mappers/blob/master/3_Quality_Preparations_Inital_FastQC.md).

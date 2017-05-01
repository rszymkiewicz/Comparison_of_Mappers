#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#This script must be started in the folder containing both VCF_File_Paths from both mappers. In this case both mapper VCF_File_Paths were moved to /media/sdc/

#Variables
SamtoolsFolderPath = "/media/sdc/mpileup/"
timestamp = (datetime.datetime.now())

print ("\n ##########\n VCF-compare (Bowtie2 and BWA mpileup) - Obtain B2 File Paths \n ##########\n")
print (timestamp)

os.system("mkdir vcf_comparison")
os.system("cd vcf_comparison")

with open('VCF_File_Path_B2.txt','rt') as B2_VCF_Path, open('VCF_File_Path_BWA.txt', 'rt') as BWA_VCF_Path:
	lineBWA = ""
	for lineB2 in B2_VCF_Path: #Bowtie2_readline: #!= lineBWA in BWA_VCF_Path:
			x = 0
			print lineB2
			B2_VCF = lineB2.replace("\n", "")
			print(B2_VCF)
			#Obtain Sample Name
			lineB2_array = B2_VCF.split("_")
			sample_Full_B2 = lineB2_array[1]
			print(sample_Full_B2)
			lineB2_array2 = sample_Full_B2.split("/")
			print(lineB2_array2)
			sample_B2 = lineB2_array2[0]
			print(sample_B2)
			print (".................")
			print ("Current Sample : " + sample_B2)
			print(B2_VCF)
			#exit()
			while x == 0:
				#read in second file
				BWA_readline = BWA_VCF_Path.readline()
				BWA_VCF = BWA_readline.replace("\n", "")
				print(BWA_VCF)
				#Obtain Sample Name
				lineBWA_array = BWA_VCF.split("_")
				sample_Full_BWA = lineBWA_array[1]
				print(sample_Full_BWA)
				lineBWA_array2 = sample_Full_BWA.split("/")
				print(lineBWA_array2)
				sample_BWA = lineBWA_array2[0]
				print(sample_BWA)
				#exit()
				if sample_B2 == sample_BWA:
					#Tabix
					print("\nTabix:")
					print("tabix -p vcf " + B2_VCF)
					os.system("tabix -p vcf " + B2_VCF)
					print("tabix -p vcf " + BWA_VCF)
					os.system("tabix -p vcf " + BWA_VCF)
					
					#vcf-compare tool in vcftools which compares 2 vcf files and outputs a Venn-diagram
					print("\nvcf-compare:")
					#Print vcf-compare code to screen as progress indication for user
					print ("vcf-compare " + B2_VCF + " " + BWA_VCF + " > /media/sdc/vcf_comparison/" + sample_B2 + ".compared")
					os.system("vcf-compare " + B2_VCF + " " + BWA_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + ".compared")
					
					#bcftools stats and bedtools plot
					print("BEDTOOLS STATS AND PLOT")
					#Obtain statistical information for Bowtie2
					bcftoolsStats = ("/home/rachelle/bcftools/bcftools stats --threads 11 -F " + referenceFasta + " -s - " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + B2_VCF " > " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + "_B2.vcf.gz.stats")
					#Print bcftools stats code to screen as progress indication for user
					print(bcftoolsStats)
					os.system(bcftoolsStats)
					plot = ("/home/rachelle/bcftools/plot-vcfstats -p " + str(Sample) + " " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + "_B2.vcf.gz.stats")
					#Print bcftools plot code to screen as progress indication for user
					print(plot)
					os.system(plot)
					#Obtain statistical information for BWA_MEM
					bcftoolsStats = ("/home/rachelle/bcftools/bcftools stats --threads 11 -F " + referenceFasta + " -s - " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + BWA_VCF " > " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + "_BWA.vcf.gz.stats")
					#Print bcftools stats code to screen as progress indication for user
					print(bcftoolsStats)
					os.system(bcftoolsStats)
					plot = ("/home/rachelle/bcftools/plot-vcfstats -p " + str(Sample) + " " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + "_BWA.vcf.gz.stats")
					#Print bcftools plot code to screen as progress indication for user
					print(plot)
					os.system(plot)
						
					#bedtools intersect tool in bedtools which identifies overlaping genomic features across 2 sets of data. In this output (defaults not -wa -wb) the output pertains to only the regions which overlap between both data sets. http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
					print("\nbedtools intersect:")
					#Print bedtools intersect code to screen as progress indication for user
					print ("bedtools intersect -a " + B2_VCF + " -b " + BWA_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + "_common_B2_BWA.vcf")
					os.system("bedtools intersect -a " + B2_VCF + " -b " + BWA_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + "_common_B2_BWA.vcf")
					
					#bedtools subtract tools in bedtools identifies features in B which overlap with A and removes the overlap reporting the remaining portion of A only.   Thus output is unique to the first parameter! http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html
					print("\nbedtools subtract:")
					#Print bedtools subtract code to screen as progress indication for user (unique Bowtie2)
					print ("bedtools subtract -a " + B2_VCF + " -b " + BWA_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + "_unique_B2.vcf")
					os.system("bedtools subtract -a " + B2_VCF + " -b " + BWA_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + "_unique_B2.vcf")
					#Print bedtools subtract code to screen as progress indication for user (unique BWA-MEM)
					print ("bedtools subtract -a " + BWA_VCF + " -b " + B2_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + "_unique_BWA.vcf")
					os.system("bedtools subtract -a " + BWA_VCF + " -b " + B2_VCF + " > /media/sdc/vcf_comparison/Sample_" + sample_B2 + "/" + sample_B2 + "_unique_BWA.vcf")
					print (".................")
					x = 1
Path_Location.close()
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
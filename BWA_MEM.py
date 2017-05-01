#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/Trimmomatic/"
timestamp = (datetime.datetime.now())
BWAFolder = "/media/sdc/BWA/Sample_"

print ("\n ##########\n BWA Alignment \n ##########\n")
print (timestamp)

#Make directory to contain all alignment files of BWA
os.system("mkdir /media/sdc/BWA/")
#Change directory to Bowtie2
os.chdir("/media/sdc/BWA/")

#BWA Library
Library = ("bwa index /media/rachelle/Seagate\ Expansion\ Drive/PES_Files/P749_PES_Sample/PES_P749_PacBio_edited.fasta")
#Print Library Build Progress to Screen
print ("Library Build: " + Library)
os.system(Library)
#Print Library Build to Screen 
print ("........... Library Complete ...........")

for root, dirs, files in os.walk(Root_Folder):
	print ("..........")
	print ("Current Directory : " + root) #root directory we are currently in
	Forward_paired = ""
	Reverse_paired = ""
	for file in files: #loop through all files within this current folder
		#Get Sample Name
		SampleName = file.split('_')
		Sample = SampleName[0]
		#Find Forward_Paired Trimmed (R1 Lane1 and R1 Lane2) files
		if file.endswith('_forward_paired.fq'):
			Forward_paired = file
		#Find Reverse_Paired Trimmed (R2 Lane1 and R2 Lane2) files
		elif file.endswith('_reverse_paired.fq'):
			Reverse_paired = file
		else:
			continue
	if Forward_paired != "" and Reverse_paired != "":
		#Make separate directory for the sample
		os.system("mkdir /media/sdc/BWA/Sample_" + str(Sample))
		#Change directory to the sample
		os.chdir("/media/sdc/BWA/Sample_" + str(Sample))
		#BWA-MEM Alignment
		BWA_align = ("bwa mem /media/rachelle/Seagate\ Expansion\ Drive/PES_Files/P749_PES_Sample/PES_P749_PacBio_edited.fasta " + (os.path.join(root, Forward_paired)) + " " + (os.path.join(root, Reverse_paired)) + " > " + BWAFolder + str(Sample) + "/" + str(Sample) + "_aln_PE.sam")
		#Print BWA-MEM alignment code to screen as progress indication for user
		print (BWA_align)
		print(os.system(BWA_align))
		#Conversion of Mapping Output File
		#Convert sam to bam
		samBam = "samtools view -S -h -F 4 -b " + BWAFolder + str(Sample) + "/" + str(Sample) + "_aln_PE.sam > " + BWAFolder + str(Sample) + "/" + str(Sample) + "_aln_PE.bam"
		#Print sam to bam code to screen as progress indication for user
		print (samBam)
		os.system(samBam)
		#Convert bam to sorted bam
		sortedBam = "samtools sort " + BWAFolder + str(Sample) + "/" + str(Sample) + "_aln_PE.bam " + BWAFolder + str(Sample) + "/" + str(Sample) + "_aln_PE.sort"
		#Print bam to sorted bam code to screen as progress indication for user
		print (sortedBam)
		os.system(sortedBam)
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
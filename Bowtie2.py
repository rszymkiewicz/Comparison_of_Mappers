#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/Trimmomatic/"
timestamp = (datetime.datetime.now())
Bowtie2Folder = "/media/sdc/Bowtie2/Sample_"

print ("\n ##########\n Bowtie2 Alignment \n ##########\n")
print (timestamp)

#Make directory to contain all alignment files of Bowtie2
os.system("mkdir /media/sdc/Bowtie2/")
#Change directory to Bowtie2
os.chdir("/media/sdc/Bowtie2/")

#Bowtie2 Library
Library = ("bowtie2-build /media/rachelle/Seagate\ Expansion\ Drive/PES_Files/P749_PES_Sample/PES_P749_PacBio_edited.fasta PacBioP749")
Library_Name = "PacBioP749"
print ("Library Build: " + Library)
os.system(Library)

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
		os.system("mkdir /media/sdc/Bowtie2/Sample_" + str(Sample))
		#Change directory to the sample
		os.chdir("/media/sdc/Bowtie2/Sample_" + str(Sample))
		#Bowtie2 Alignment
		bowtie2_align = ("bowtie2 -t -p 10 -x /media/sdc/Bowtie2/" + Library_Name + " -1 " + (os.path.join(root, Forward_paired)) + " -2 " + (os.path.join(root, Reverse_paired)) + " -S " + Bowtie2Folder + str(Sample) + "/" + str(Sample) + "_Reads_to_Contig.sam --no-unal --no-mixed")
		#Print Bowtie2 alignment code to screen as progress indication for user
		print (bowtie2_align)
		#Perform bowtie2 mapping
		os.system(bowtie2_align)
		#Conversion of Mapping Output File
		#Convert sam to bam
		samBam = "samtools view -S -h -F 4 -b " + Bowtie2Folder + str(Sample) + "/" + str(Sample) + "_Reads_to_Contig.sam > " + Bowtie2Folder + str(Sample) + "/" + str(Sample) + "_Reads_to_Contig.bam"
		#Print sam to bam code to screen as progress indication for user
		print (samBam)
		os.system(samBam)
		#Convert bam to sorted bam
		sortedBam = "samtools sort " + Bowtie2Folder + str(Sample) + "/" + str(Sample) + "_Reads_to_Contig.bam " + Bowtie2Folder + str(Sample) + "/" + str(Sample) + "_Reads_to_Contig.sort"
		#Print bam to sorted bam code to screen as progress indication for user
		print (sortedBam)
		os.system(sortedBam)
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
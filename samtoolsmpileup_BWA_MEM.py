#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/BWA/"
referenceFasta = "/media/rachelle/Seagate\ Expansion\ Drive/PES_Files/P749_PES_Sample/PES_P749_PacBio_edited.fasta"
SamtoolsFolderPath = "/media/sdc/mpileupBWA/"
sortedBam = ""
timestamp = (datetime.datetime.now())

print ("\n ##########\n samtools mpileup BWA \n ##########\n")
print (timestamp)

#Make Parent Directory
os.system("mkdir /media/sdc/mpileupBWA/")

for root, dirs, files in os.walk(Root_Folder):
	print ("..........")
	print ("Current Directory : " + root) #root directory we are currently in
	Forward_paired = ""
	Reverse_paired = ""

	for file in files: #loop through all files within this current folder
		#Get Sample Name
		SampleName = file.split('_')
		Sample = SampleName[0]
				
		#Find the sorted.bam file from the BWA-MEM mapping
		if file.endswith('.sort.bam'):
			sortedBam = file	
		else:
			continue
	if sortedBam != "":
		#Make separate directory for the sample
		#print ("mkdir " + SamtoolsFolderPath + "Sample_" + str(Sample))
		os.system("mkdir " + SamtoolsFolderPath + "Sample_" + str(Sample))
		#Change directory to the sample
		os.chdir(SamtoolsFolderPath + "Sample_" + str(Sample))
		#Sample Index
		print("SAMPLE INDEX")
		IndexBuild = ("samtools faidx " + referenceFasta)
		#Print samtools index code to screen as progress indication for user
		print(IndexBuild)
		os.system(IndexBuild)
		print("..........")
		print("..........")
		#mpileup
		print("MPILEUP")
		#Create mpileup VCF file per sample
		mpileup = ("/home/rachelle/samtools/samtools mpileup -vf " + referenceFasta + " " + os.path.join(root, sortedBam) + " | /home/rachelle/bcftools/bcftools call -vmO z -o " + SamtoolsFolderPath + "Sample_" + str(Sample) + "/" + str(Sample) + "_BWA.vcf.gz") 
		#Print samtools mpileup code to screen as progress indication for user
		print(mpileup)
		os.system(mpileup)
		print("..........")
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
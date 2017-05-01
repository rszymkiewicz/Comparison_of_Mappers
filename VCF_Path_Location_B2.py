#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/mpileup/"
referenceFasta = "/media/rachelle/Seagate\ Expansion\ Drive/PES_Files/P749_PES_Sample/PES_P749_PacBio_edited.fasta"
SamtoolsFolderPath = "/media/sdc/mpileup/"
sortedBam = ""
timestamp = (datetime.datetime.now())

print ("\n ##########\n VCF-Compare Preparation (Bowtie2 and BWA mpileup) - Obtain B2 File Paths \n ##########\n")
print (timestamp)
#find /media/sdc/mpileup/ -type f -name "*_B2.vcf.gz" > VCF_File_Path_B2.txt


with open('VCF_File_Path_B2.txt','w') as Path_Location:
		
	for root, dirs, files in os.walk(Root_Folder):
		print ("..........")
		print ("Current Directory : " + root) #root directory we are currently in
		
		for file in files: #loop through all files within this current folder
			#Get Sample Name
			SampleName = file.split('_')
			Sample = SampleName[0]
					
			#Find VCF Output file from Bowtie2 samtools mpileup 
			if file.endswith('.vcf.gz'):
				vcf = file
				#write the file path of the Bowtie2 VCF file to a new file
				Path_Location.write(SamtoolsFolderPath + vcf + "\n")
			else:
				continue	
Path_Location.close()
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
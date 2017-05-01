#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/mpileupBWA/"
referenceFasta = "/media/rachelle/Seagate\ Expansion\ Drive/PES_Files/P749_PES_Sample/PES_P749_PacBio_edited.fasta"
SamtoolsFolderPath = "/media/sdc/mpileupBWA/"
sortedBam = ""
timestamp = (datetime.datetime.now())

print ("\n ##########\n VCF-compare (Bowtie2 and BWA mpileup) - Obtain BWA File Paths \n ##########\n")
print (timestamp)
#find /media/sdc/mpileupBWA/ -type f -name "*_BWA.vcf.gz" > VCF_File_Path_BWA.txt


with open('VCF_File_Path_BWA.txt','w') as Path_Location:
		
	for root, dirs, files in os.walk(Root_Folder):
		print ("..........")
		print ("Current Directory : " + root) #root directory we are currently in
		Forward_paired = ""
		Reverse_paired = ""

		for file in files: #loop through all files within this current folder
			#Get Sample Name
			SampleName = file.split('_')
			Sample = SampleName[0]
					
			#Find VCF Output file from BWA-MEM samtools mpileup
			if file.endswith('.vcf.gz'):
				vcf = file
				#write the file path of the BWA-MEM VCF file to a new file
				Path_Location.write(SamtoolsFolderPath + vcf + "\n")
			else:
				continue	
Path_Location.close()
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
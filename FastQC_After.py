#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/Trimmomatic/"
timestamp = (datetime.datetime.now())

print ("\n ##########\n Quality Check 2 [After]: FastQC Forward_Paired and Reverse_Paired \n ##########\n")
print (timestamp)

#Make directory to contain all FastQC After Trimming results
os.system("mkdir /media/sdc/fastqc_2")

for root, dirs, files in os.walk(Root_Folder):
	print ("..........")
	print ("Current Directory : " + root) #root directory we are currently in
	for file in files: #loop through all files within this current folder
		#Get Sample Name
		SampleName = file.split('_')
		Sample = SampleName[0]
		#Make separate directory for the sample
		os.system("mkdir /media/sdc/fastqc_2/Sample_" + str(Sample))
		#Check to see Forward_Paired Trimmed file
		if file.endswith('_forward_paired.fq.gz'):
			R1 = file
			output_R1 = ("fastqc " + (os.path.join(root, R1)) + " -o /media/sdc/fastqc_2/Sample_" + str(Sample))
			#Perform FastQC on the Forward_Paired sample
			os.system(output_R1)
			#Print progress to the screen that Forward_Paired FastQC was complete
			print ('\t' + 'Forward paired Fastqc Complete')
		#Check to see Reverse_Paired Trimmed file	
		if file.endswith('_reverse_paired.fq.gz'):
			R2 = file
			output_R2 = ("fastqc " + (os.path.join(root, R2)) + " -o /media/sdc/fastqc_2/Sample_" + str(Sample))
			#Perform FastQC on the Reverse_Paired sample
			os.system(output_R2)
			#Print progress to the screen that Reverse_Paired FastQC was complete
			print ('\t' + 'Reverse paired Fastqc Complete')
		else: 
			continue
	#Print dividing line between Samples
	print ('..........')
		
timestamp = (datetime.datetime.now())		
print ("Script is complete")
print (timestamp)
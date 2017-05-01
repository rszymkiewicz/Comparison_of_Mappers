#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
Root_Folder = "/media/sdc/raw_seq/"
timestamp = (datetime.datetime.now())

print ("\n ##########\n Quality Check [Before]: FastQC R1 and R2 \n ##########\n")
print (timestamp)

for root, dirs, files in os.walk(Root_Folder):
	print ("..........")
	print ("Current Directory : " + root) #root directory we are currently in
	for file in files: #loop through all files within this current folder
		#Get Sample Name
		SampleName = file.split('_')
		Sample = SampleName[0]
		#Make separate directory for the sample
		os.system("mkdir /media/sdc/fastqc_1/Sample_" + str(Sample))
		#Check to see R1 concatenated file
		if file.endswith('_L00X_R1_00X.fastq'):
			R1 = file
			output_R1 = ("fastqc " + (os.path.join(root, R1)) + " -o /media/sdc/fastqc_1/Sample_" + str(Sample))
			#Perform FastQC on the R1 sample
			os.system(output_R1)
			#Print progress to the screen that R1 FastQC was complete
			print ('\t' + 'R1 Fastqc Complete')
		#Check to see R2 concatenated file
		if file.endswith('_L00X_R2_00X.fastq'):
			R2 = file
			#Perform FastQC on the R2 sample
			output_R2 = ("fastqc " + (os.path.join(root, R2)) + " -o /media/sdc/fastqc_1/Sample_" + str(Sample))
			os.system(output_R2)
			#Print progress to the screen that R2 FastQC was complete
			print ('\t' + 'R2 Fastqc Complete')
		else: 
			continue
	#Print dividing line between Samples		
	print ('..........')
		
timestamp = (datetime.datetime.now())		
print ("Script is complete")
print (timestamp)	
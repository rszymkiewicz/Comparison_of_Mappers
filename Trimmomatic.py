#!/usr/bin/python
#Rachelle Szymkiewicz

#Libraries
import os
import datetime

#Variables
R1 = ""
R2 = ""
Root_Folder = "/media/sdc/raw_seq"
timestamp = (datetime.datetime.now())

print ("\n ##########\n Quality Check [Trimming of Reads]: Trimmomatic R1 and R2 \n ##########\n")
print (timestamp)

#Make new directory for Trimming
os.system("mkdir /media/sdc/Trimmomatic")

for root, dirs, files in os.walk(Root_Folder):
	print ("..........")
	print ("Current Directory : " + root) #root directory we are currently in		
	for file in files: #loop through all files within this current folder
		#Get Sample Name
		SampleName = file.split('_')
		Sample = SampleName[0]
		#Make separate directory for the sample
		os.system("mkdir /media/sdc/Trimmomatic/Sample_" + str(Sample))
		#Check to see R1 concatenated file
		if file.endswith('_L00X_R1_00X.fastq'):
			R1 = file
		#Check to see R2 concatenated file
		if file.endswith('_L00X_R2_00X.fastq'):
			R2 = file
		else: 
			continue
	if R1 != "" and R2 != "":
		#Change the directory to the sample directory listed
		os.chdir("/media/sdc/Trimmomatic/Sample_" + str(Sample))		
		trimmomatic_cmd = ("java -jar /home/rachelle/bin/trimmomatic.jar PE -threads 12 -trimlog " + Sample + " " + (os.path.join(root, R1)) + " " + (os.path.join(root, R2)) + " " + Sample + "_forward_paired.fq.gz " + Sample +"_forward_unpaired.fq.gz " + Sample + "_reverse_paired.fq.gz " + Sample +"_reverse_unpaired.fq.gz " + "ILLUMINACLIP:/home/rachelle/bin/trimmomatic-master/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
		#Perform Trimmomatic on R1 and R2 of sample
		#print (trimmomatic_cmd)
		os.system(trimmomatic_cmd)
		#Print progress to the screen that trimming of sample was complete
		print ('Trim Complete')
	else: 
		#If R1 and R2 files were not able to be found for the sample, no trimming would occur and an output message would appear on screen
		print ('No trimming')
timestamp = (datetime.datetime.now())		
print ('Script is complete')
print (timestamp)
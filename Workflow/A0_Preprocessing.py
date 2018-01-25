#!/usr/bin/env python

### 25/01/2018
### In house RNAseq workflow for the de Bruijn lab

#####################

import argparse, os, csv, glob, sys

outf = open('A0_Progress.txt', 'w')

"""
Argument Parsing
"""

parser = argparse.ArgumentParser(
	description='WIP basic FastQ QC and pre-processing pipeline. Trigger pipeline in the current directory you want the output to be generated.'
)

parser.add_argument(
	'-f',
	help = 'Required. Directory containing FastQ files',
	required = True
)
parser.add_argument(
	'-a',
	choices=['nextera', 'illumina', 'small_rna'],
	help = 'Specify adaptor used during library prep (nextera, illumina, small_rna). Default is auto detect, when not specified'
)
parser.add_argument(
	'-s',
	help = 'Required. Specify the Samples.txt file to be used. This is a tab deliminated file that specifies files, samples, reps, lanes and pairing.',
	required = True
)
parser.add_argument(
	'-namechange',
	help = 'True if you want file names to be changed based on sample annotations. Niche use, for when files are already the appropriate name.'
)
parser.add_argument(
	'-trim',
	help = "Required. True or False. Specify whether files should be trimmed. Trimming is recommended, unless the trimming output is concerning.",
	required = True
)
parser.add_argument(
	'-maprun',
	help = "Required. True or False. Specify whether mapping should be run in pipeline. Set this to False if you have too many files, or want to verify the commands beforehand.",
	required = True
)
parser.add_argument(
        '--ERCC',
        help = 'Flag for ERCC. If present, will map reads to mm9 + ERCC',
        action = 'store_true'
)
parser.add_argument(
	'--FC',
	help = 'Flag for Feature Counts. If present, will use subread tool to count reads mapping to genetic exons. If skipped, will still produce .sh files to generate Feature Counts table independently.',
        action = 'store_true'
)
parser.add_argument(
	'--PCR_Dedup',
	help = 'Flag for removing PCR duplicates after mapping. If present, will use samtools to remove any PCR duplicates found in the .BAM files.',
        action = 'store_true'
)
parser.add_argument(
	'--BW',
	help = 'Flag for processing BigWigs. If present, will create bigwig files for displaying on UCSC.',
        action = 'store_true'
)

args = parser.parse_args()

"""
Output preamble
"""

### This should be more verbose, and should explain all parameters in use. ###
outf.write("Mapping pipeline\n\nParameters:\n")
outf.write("FastQ location:\t" + args.f)
outf.write("Samples.txt:\t" + args.s)
outf.write("Trimming:\t" + args.trim)
outf.write("File name change?:\t" + args.namechange)
outf.write("Run mapping?:\t" + args.maprun)

"""
Sample setup
   - Requires Sample.txt to be set up correctly!
"""

outf.write("--------------------\nA0_Preprocessing pipeline\n--------------------\n\n")

# Read Sample sheet

filenames=[]
samples=[]
reps=[]
lanes=[]
pairs=[]

with open(args.s,'r') as f:
    next(f) # skip headings
    reader=csv.reader(f,delimiter='\t')
    for A,B,C,D,E in reader:
        filenames.append(A)
        samples.append(B) 
	reps.append(C)
	lanes.append(D)
	pairs.append(E)

# Determine whether data is paired

if len(set(pairs)) == 2:
	paired = True
elif len(set(pairs)) == 1:
	paired = False
else:
	outf.close("ERROR: Invalid number of pair values, cannot determine pairing")
	sys.exit()

# Summarise data and print

samplesSet = set(samples)
repsNum = len(set(reps))
lanesNum = len(set(lanes))
fileNum = len(filenames)
	
outf.write("Description of the data:\n\n\tNumber of files: " + str(fileNum) + "\n\tNumber of reps: " + str(repsNum) + "\n\tNumber of lanes: " + str(lanesNum) + "\n\tSample groups: " + str(samplesSet) + "\n\n")

"""
Error handling!
"""

# Variable lengths consistent?
if (len(filenames) != len(samples) or len(filenames) != len(reps) or len(filenames) != len(lanes) or len(filenames) != len(pairs)):
	outf.write("ERROR: Number of items not consistent in setup file")
	sys.exit()
	
# Checking if pair values 1 or 2
for i in pairs:
	if (i != '1' and i != '2'):
		outf.write("ERROR: Pair values must only contain 1 or 2")
		sys.exit()

# Equal number of paired files?
if paired == True:
	if pairs.count('1') != pairs.count('2'):
		outf.write("ERROR: Number of pairs do not match.")
		sys.exit()
		
# Are all filenames unique?
if len(filenames) != len(set(filenames)):
	outf.write("ERROR: Not all filenames are unique")
	sys.exit()

"""
Create tmp files and rename
"""

# New file names!

### Currently handles .gz compression and .sanfastq.gz, and .fastq/sanfastq
### Any other file types to be used must be specified here, or a cleaner solution found.

newFilenames = []
extension = ""

for f in range(fileNum):
	if filenames[f][-3:] == ".gz":
		if filenames[f][-11:] == "sanfastq.gz":
			extension = ".sanfastq.gz"
		else:
			extension = ".fastq.gz"
	else:
		if filenames[f][-8:] == "sanfastq":
			extension = ".sanfastq"
		else:
			extension = ".fastq"
	fn = samples[f] + "_" + reps[f] + "_" + lanes[f] + "_" + pairs[f] + extension
	newFilenames.append(fn)


# Rename old files

outf.write("Creating tmp Fastq files and renaming to more appropriate labels\n")

os.system("mkdir -p ./tmpFastqFiles")
os.system("cp " + args.f +"/*fastq* ./tmpFastqFiles")
os.system("cp " + args.f +"/*.fq* ./tmpFastqFiles")
 
if args.namechange == 'True':
	for fq in glob.iglob("./tmpFastqFiles/*"):
		ind = filenames.index(os.path.basename(fq))
		newName = newFilenames[ind]
		os.system("mv " + fq + " ./tmpFastqFiles/" + newName)

fqDir = "./tmpFastqFiles"

""" 
FastQC - first run
   Establish pre-trimming quality
"""
outf.write("Running FastQC...\n")

if args.trim == "True":
	os.system("mkdir -p ./FastQC_PreTrim")
	CmdFastQC = "fastqc -o ./FastQC_PreTrim -t 2 " + fqDir + "/*fastq*"
else:
	os.system("mkdir -p ./FastQC")
	CmdFastQC = "fastqc -o ./FastQC -t 2 " + fqDir + "/*fastq*"	
os.system(CmdFastQC)

outf.write("FastQC complete!\n\n")


"""
TrimGalore
   Trim for quality and adaptors
   Code split based off presence of argument -a and pairing
"""

# Build trim command

### Could do with optimising trim_galore parameters.
### Could also attempt repeated mapping runs, with progressively aggressive trimming and taking the optimimum mapping rate - Very work intensive!

if args.trim == "True":
	outf.write("Trimming adaptors...\n")
	
	os.system("mkdir -p ./trim")
	
	if args.a:
		if paired == True:
			outf.write("\tRunning TrimGalore for " + args.a + " adaptors, with paired end reads\n")
			CmdTrim = "trim_galore --paired --" + args.a + " --stringency 3 --output_dir ./trim " + fqDir + "/*fastq*"
		else:
			outf.write("\tRunning TrimGalore for " + args.a + " adaptors, with single end reads\n")
			CmdTrim = "trim_galore --" + args.a + " --stringency 3 --output_dir ./trim " + fqDir + "/*fastq*"		
	else:
		if paired == True:
			outf.write("\tRunning TrimGalore with adaptor auto-detection. with paired end reads\n")
			CmdTrim = "trim_galore --paired --stringency 3 --output_dir ./trim " + fqDir + "/*fastq*"
		else:
			outf.write("\tRunning TrimGalore with adaptor auto-detection. with single end reads\n")
			CmdTrim = "trim_galore --stringency 3 --output_dir ./trim " + fqDir + "/*fastq*"
			
	outf.write("\tCommand:\t" + CmdTrim + "\n")
	os.system(CmdTrim)
	os.system('mkdir -p ./trim/trimOutput')
	os.system('mv ./trim/*.txt ./trim/trimOutput')
	
	outf.write("\tAdaptors trimmed, Fastq files stored in ./trim, Trim stats stored in ./trim/trimOutput\n")
	
	# Call perl script for summarising Trimming statistics
	
	outf.write("\tExtracting summary of trim statistics to outputTrimStats.txt\n")
	os.system("perl /t1-data/user/jharman/pipeline_v2/A1_Gather_Trim_statistics.pl ./trim/trimOutput/*.txt")
	
	# Call R script for visualising Trimming statistics
	
	outf.write("\tVisualising FastQ and Trimming statistics\n")
	os.system("Rscript /t1-data/user/jharman/pipeline_v2/A2_TrimPlot.R")
	
	fqDir = "./trim"  # Change fqDir to trim only if trimming occured!
else:
	outf.write("Trimming skipped\n\n")
	fqDir = "./tmpFastqFiles"	


"""
FastQC - second run
   Post-trimming quality
"""
if args.trim == "True":
	outf.write("Running FastQC on trimmed FastQ...\n")

	os.system("module load fastqc/0.11.4")
	os.system("mkdir -p ./FastQC_PostTrim")

	CmdFastQC = "fastqc -o ./FastQC_PostTrim -t 2 ./trim/*fq*"
	os.system(CmdFastQC)

	outf.write("Trimmed FastQC complete!\n\n")

outf.write("Pre-processing pipeline complete!\n\n")

outf.write("--------------------\nB0_Mapping pipeline\n--------------------\n\n")

"""
Rebuild new trimmed file names!
"""

# Required, as trim galore changes file names, and changes them based on whether the data is paired or single end.

trimFilenames = []

for f in range(fileNum):
	if args.trim == "True":
        	if args.namechange == "True":
			if len(set(pairs)) == 2:
				if pairs[f] == '1':
                                	fn = newFilenames[f].split(os.extsep)[0]+'_val_1.fq.gz'
                                	trimFilenames.append(fn)
                        	elif pairs[f] == '2':
                                	fn = newFilenames[f].split(os.extsep)[0]+'_val_2.fq.gz'
                                	trimFilenames.append(fn)
                        	else:
                                	outf.write("ERROR: Pairs must be 1 or 2. Something has gone terribly wrong.")
                                	sys.exit()
                	elif len(set(pairs)) == 1:
                        	fn = newFilenames[f].split(os.extsep)[0]+'_trimmed.fq.gz'
                        	trimFilenames.append(fn)
        	else:
                	if len(set(pairs)) == 2:
                        	if pairs[f] == '1':
                                	fn = filenames[f].split(os.extsep)[0]+'_val_1.fq.gz'
                                	trimFilenames.append(fn)
                        	elif pairs[f] == '2':
                                	fn = filenames[f].split(os.extsep)[0]+'_val_2.fq.gz'
                                	trimFilenames.append(fn)
                        	else:
                                	outf.write("ERROR: Pairs must be 1 or 2. Something has gone terribly wrong.")
                                	sys.exit()
                	elif len(set(pairs)) == 1:
                        	fn = filenames[f].split(os.extsep)[0]+'_trimmed.fq.gz'
                        	trimFilenames.append(fn)
	else:
		fn = newFilenames[f]
		trimFilenames.append(fn)
			
"""
Constructing STAR mapping command
"""

# This is a complicated loop, but it works to construc the Star command! Leave as is, unless you want to optimise code.

sampleCount = 0
repCount = 0
pairCount = 0
laneCount = 0
os.system("mkdir -p ./StarCommands")
os.system("mkdir -p ./StarBAM")
os.system("module load rna-star")
outf.write("Seeding files for mapping...\n\n")

# These first two loops are the seeding information to match the file being mapped. - i.e, the sample name and rep number.

for s in range(len(set(samples))):
        for r in range(len(set(reps))):
                cmd = []  # List of files for each sample mapping
                suffix = ""  # Mapping suffix for each sample mapping
                starName= ""  # File name for each sample mapping
		
                # These two loops specify the pair and lane seeding required for mapping the individual file (e.g, 3 lanes, each with 2 pairs)
                for p in range(len(set(pairs))):
                        for l in range(len(set(lanes))):
                                
                                # This next code block is to loop through every file until we come across the one file that matches everything as seeded above!
				# i.e, loops over every file to find what matches the sample name, rep number, lane number, and pairing. Then generates the sample suffix and file names.
				# note that the order of these loops is very deliberate - the structure of the command is as follows: lane1pair1,lane2pair1,lane3pair1 lane1pair2,lane2pair2,lane3pair2
                                for i in range(len(newFilenames)):
                                        if samples[i] == list(set(samples))[s]:
                                                sampleCount +=1
                                                if reps[i] == list(set(reps))[r]:
                                                        repCount +=1
                                                        if pairs[i] == list(set(pairs))[p]:
                                                                pairCount +=1
                                                                if lanes[i] == list(set(lanes))[l]:
                                                                        # Once a match is acquired, append trim file name, and create file suffix
                                                                        laneCount +=1
                                                                        cmd.append(trimFilenames[i])
                                                                        suffix = samples[i] + "_" + reps[i] + "_"
                                                                        starName = "./StarCommands/" + suffix + "Mapping.sh"
                
                # Now that we have a list of filenames in the correct order, we can string them together to create a command!
                # We also create the .sh file here
                if suffix != "":
                        Command = ""
                        if paired == True:
                                for comm in range(len(cmd)):
                                        x = len(cmd)/2-1  # this is a quick calc to determine the half way file that changes the concatenation from "," to " "
					if args.trim == "True":
						fileLoc = "./trim/" + cmd[comm]
                                        	Command += fileLoc
					else:
						fileLoc = "./tmpFastqFiles/" + cmd[comm]
						Command += fileLoc
                                        if comm == len(cmd)-1:
                                                break
                                        if comm < x:
                                                Command += ","
                                        elif comm == x:
                                                Command += " "  # This is the half way concatenation. Changeover here from pair 1 to pair 2.
                                        elif comm > x:
                                                Command += ","
                        if paired == False:
                                for comm in range(len(cmd)):
					if args.trim == "True":
						fileLoc = "./trim/" + cmd[comm]
                                        	Command += fileLoc
					else:
						fileLoc = "./tmpFastqFiles/" + cmd[comm]
						Command += fileLoc
                                        if comm == len(cmd)-1:
                                                break
                                        if comm < len(cmd)-1:
                                                Command += ","  
                        if args.ERCC == True:
                                StarSH = open(starName, 'w')
                                StarSH.write("#!/bin/sh\n\n#$ -cwd\n\n#$ -q batchq\n\n#$ -M " + args.user + "\n#$ -m eas\n\nmodule load rna-star\n\n")
                                StarSH.write("STAR --genomeDir  /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR_with_ERCC --readFilesIn " + Command + " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./StarBAM/" + suffix + " --runThreadN 4")
                                StarSH.close()
				if args.maprun == "True":
					os.system("STAR --genomeDir  /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR_with_ERCC --readFilesIn " + Command + " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./StarBAM/" + suffix + " --runThreadN 4")
                        		outf.write("Mapping complete for " + suffix + "\n\n")
			else:
                                StarSH = open(("./StarCommands/" + suffix + "Mapping.sh"), 'w')
                                StarSH.write("#!/bin/sh\n\n#$ -cwd\n\n#$ -q batchq\n\n#$ -M " + args.user + "\n#$ -m eas\n\nmodule load rna-star\n\n")
                                StarSH.write("STAR --genomeDir  /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR --readFilesIn " + Command + " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./StarBAM/" + suffix + " --runThreadN 4")
				StarSH.close()
				if args.maprun == "True":
					os.system("STAR --genomeDir  /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR --readFilesIn " + Command + " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./StarBAM/" + suffix + " --runThreadN 4")
					outf.write("Mapping complete for " + suffix + "\n\n")
"""
Post mapping cleanup & statistics
"""
os.system("mkdir -p ./StarBAM/Log")
os.system("mv ./StarBAM/*Log.out ./StarBAM/Log")
os.system("mv ./StarBAM/*Log.final.out ./StarBAM/Log")
os.system("mv ./StarBAM/*Log.progress.out ./StarBAM/Log")
os.system("mv ./StarBAM/*SJ.out.tab ./StarBAM/Log")

# Perl and R scripts for summarising and visualising mapping
os.system("perl /t1-data/user/jharman/pipeline_v2/B1_Gather_Map_statistics.pl ./StarBAM/Log/*.final.out")
os.system("Rscript /t1-data/user/jharman/pipeline_v2/B2_MapPlot.R")


"""
Removing PCR duplicates
"""
### Could switch this with Picard if required - Picard is more intensive, but also more accurate

if args.PCR_Dedup == True:
	if paired == True:
		cmd = "samtools rmdup "
	else:
		cmd = "samtools rmdup -s "
	if args.maprun == "True":
		os.system("mkdir -p ./StarBAM/filtered")
		fqPost = ""
		for fq in glob.iglob("./StarBAM/*.bam"):
			fqPost = "./StarBAM/filtered/" + os.path.basename(fq) + ".filtered.bam"
			
			cmd2 = cmd + fq + " " + fqPost
			
			os.system(cmd2)
		outf.write("PCR duplicates removed from BAM files\n\n")
	

"""
FeatureCounts
"""

outf.write("Generating featureCounts.sh - submit this after pipeline!\n\n")

if args.PCR_Dedup == True:
	fcFiles = ' '.join(glob.glob('./StarBAM/filtered/*.bam'))
else:
	fcFiles = ' '.join(glob.glob('./StarBAM/*.bam'))

os.system("mkdir -p ./featureCounts")

if args.ERCC == True:
        if paired == True:
                FC = open("featureCounts.sh", 'w')
                FC.write("#!/bin/sh\n\n#$ -cwd\n\n#$ -q batchq\n\n#$ -M " + args.user + "\n#$ -m eas\n\nmodule load subread\n\n")
                FC.write("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -p -s 0 " + fcFiles)
                FC.close()
		if args.FC == True:
			os.system("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -p -s 0 " + fcFiles)
			outf.write("Feature Counts table Generated")
        else:
                FC = open("featureCounts.sh", 'w')
                FC.write("#!/bin/sh\n\n#$ -cwd\n\n#$ -q batchq\n\n#$ -M " + args.user + "\n#$ -m eas\n\nmodule load subread\n\n")
                FC.write("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -s 0 " + fcFiles)
                FC.close()
		if args.FC == True:
			os.system("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -s 0 " + fcFiles)
			outf.write("Feature Counts table Generated")
else:
        if paired == True:
                FC = open("featureCounts.sh", 'w')
                FC.write("#!/bin/sh\n\n#$ -cwd\n\n#$ -q batchq\n\n#$ -M " + args.user + "\n#$ -m eas\n\nmodule load subread\n\n")
                FC.write("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -p -s 0 " + fcFiles)
                FC.close()
		if args.FC == True:
			os.system("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -p -s 0 " + fcFiles)
			outf.write("Feature Counts table Generated")
        else:
                FC = open("featureCounts.sh", 'w')
                FC.write("#!/bin/sh\n\n#$ -cwd\n\n#$ -q batchq\n\n#$ -M " + args.user + "\n#$ -m eas\n\nmodule load subread\n\n")
                FC.write("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -s 0 " + fcFiles)
                FC.close()
		if args.FC == True:
			os.system("featureCounts -T 5 -a /databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf -t exon -g gene_id -o ./featureCounts/features_counted.txt -s 0 " + fcFiles)
			outf.write("Feature Counts table Generated")

"""
BigWigs
"""
if args.BW == True:
	outf.write("\nBigWig files generated\n")
	os.system("bamToGBrowse.pl ./StarBAM /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa")


"""
Closing
"""
outf.write("Removing temporary files...\n")
os.system("rm ./tmpFastqFiles/*fastq*")
os.system("rmdir ./tmpFastqFiles")

outf.write("Pipeline complete!")

outf.close()


#!/bin/sh
##########################################################################
## A script template for submitting batch jobs. To submit a batch job, 
## please type
##
##    qsub myprog.sh
##
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################

## The following to run programs in the current working directory

#$ -cwd


## Specify a queue

#$ -q batchq


## The following two lines will send an email notification when the job is 
## Ended/Aborted/Suspended - Please replace "UserName" with your own username.

#$ -M jharman
#$ -m eas

module load python
module load fastqc/0.11.4
module load trim_galore/0.4.1
module load rna-star
module load ucsctools
python /t1-data/user/jharman/pipeline_v2/A0_Preprocessing.py -user jharman -f /t1-data/user/jharman/YS_Exp/RNA/fq_files/ -s /t1-data/user/jharman/YS_Exp/RNA/Samples.txt -namechange True -trim True -maprun True

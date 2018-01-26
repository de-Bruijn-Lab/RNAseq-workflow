#!/bin/sh

#$ -cwd
#$ -q batchq
#$ -M username
#$ -m eas

module load python
module load fastqc/0.11.4
module load trim_galore/0.4.1
module load rna-star
module load ucsctools

python /t1-data/user/jharman/Scripts/RNAseq/RNAseqWorkflow.py -user jharman -f /path/to/fq.files/ -s /path/to/SampleSheet.txt -namechange True -trim True -maprun True --ERCC

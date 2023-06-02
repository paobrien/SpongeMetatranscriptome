#!/bin/sh

### CONCATENATE SEQUENCES FROM THE SAME SAMPLE IN DIFFERENT RUNS/LANES ###

# use this script to submit each folder to qsub to run in parallel

for PREF in $(ls /scratch/jc254091/Metagenomics/Fastq_concatenated/); do
# Check variables
  echo "prefix is $PREF"
  echo "file path is /scratch/jc254091/Metagenomics/Fastq_concatenated/${PREF}"
  ls /scratch/jc254091/Metagenomics/Fastq_concatenated/${PREF}
# Submit to qsub to run in parallel
  qsub 01_batch_concatenate.pbs -v "DIR=/scratch/jc254091/Metagenomics/Fastq_concatenated/${PREF}/"
done


# END SCRIPT






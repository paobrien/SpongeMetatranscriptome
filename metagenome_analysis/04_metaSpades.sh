#!/bin/sh

# Use this script to submit jobs in parallel via qsub

for FILE in $(ls /scratch/jc254091/Metagenomics/SOAPnuke/To_assemble/*_1_soapnuke.fq.gz); do
  PREFIX=$(basename $FILE | cut -d "_" -f 1)
  DIR=$(dirname $FILE)
  echo "For the file $FILE, the prefix is $PREFIX"
  echo "The root path is $DIR"
  echo "The forward reads to submit are ${PREFIX}_1_soapnuke.fq.gz"
  echo "The reverse reads to submit are ${PREFIX}_2_soapnuke.fq.gz"
  qsub ~/Scripts/Metagenomics/04_metaSpades.pbs -v "FORWARD=${DIR}/${PREFIX}_1_soapnuke.fq.gz,REVERSE=${DIR}/${PREFIX}_2_soapnuke.fq.gz"
done






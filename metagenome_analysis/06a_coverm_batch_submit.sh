#!/bin/sh

##############################################
### BATCH SUBMIT SCAFFOLDS TO MAKE BAM FILES #
##############################################

# Edit files and qsub script as appropriate

# --------------------------------------------


## ALL READS AGAINST ALL SCAFFOLDS ##


#for FILE in $(ls /scratch/jc254091/Metagenomics/MetaSpades/506_spades/scaffolds.fasta); do
# Check variables
#  echo "submitting scaffold $FILE for bam alignment"

# Submit files to qsub script
#  qsub ~/Scripts/Metagenomics/06b_coverm_make_all.pbs -v "FILE=$FILE" 
#done



## FOR CARTERIOSPONGIA ##


#for FILE in $(ls /scratch/jc254091/Metagenomics/MetaSpades/57[6,7,8]_spades/scaffolds.fasta) \
#$(ls /scratch/jc254091/Metagenomics/MetaSpades/50[5,7]_spades/scaffolds.fasta); do
# Check variables
#  echo "submitting scaffold $FILE for bam alignment"

# Submit files to qsub script
#  qsub ~/Scripts/Metagenomics/06b_coverm_make_carterio.pbs -v "FILE=$FILE" 
#done



## FOR IBLUE ##


#for FILE in $(ls /scratch/jc254091/Metagenomics/MetaSpades/50[1,3,4]_spades/scaffolds.fasta) \
#$(ls /scratch/jc254091/Metagenomics/MetaSpades/516_spades/scaffolds.fasta) \
#$(ls /scratch/jc254091/Metagenomics/MetaSpades/567_spades/scaffolds.fasta); do
# Check variables
#  echo "submitting scaffold $FILE for bam alignment"

# Submit files to qsub script
#  qsub ~/Scripts/Metagenomics/06b_coverm_make_iblue.pbs -v "FILE=$FILE" 
#done


## FOR IRCINIA RAMOSA

# Note: only submitting a few samples at a time

#for FILE in $(ls /scratch/jc254091/Metagenomics/MetaSpades/50[2,6,8,9]_spades/scaffolds.fasta) \
#$(ls /scratch/jc254091/Metagenomics/MetaSpades/51[0-5]_spades/scaffolds.fasta) \
#$(ls /scratch/jc254091/Metagenomics/MetaSpades/575_spades/scaffolds.fasta); do
# Check variables
#  echo "submitting scaffold $FILE for bam alignment"

# Submit files to qsub script
#  qsub ~/Scripts/Metagenomics/06b_coverm_make_iramosa.pbs -v "FILE=$FILE" 
#done

# redo two samples that didn't finish..
for FILE in $(ls ~/Metagenomics/Assembly_and_binning/MetaSpades/515_spades/scaffolds.fasta) \
$(ls ~/Metagenomics/Assembly_and_binning/MetaSpades/575_spades/scaffolds.fasta); do
# Check variables
  echo "submitting scaffold $FILE for bam alignment"

# Submit files to qsub script
  qsub ~/Scripts/Metagenomics/06b_coverm_make_iramosa.pbs -v "FILE=$FILE" 
done




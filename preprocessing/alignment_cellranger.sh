#! /bin/bash
#$ -e ./log
#$ -o ./log
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l mem_free=150G
#$ -l h_rt=30:00:00
#$ -l scratch=200G
#$ -t 1-3:1
           
#cellranger alignment

source ~/.bashrc

#REF=/wynton/group/pollen/genomes/cellranger/GRCh38_and_rheMac10
REF=/wynton/group/pollen/genomes/cellranger/refdata-gex-GRCh38-2024-A
JOB_DIR=/wynton/scratch/jding/Cellranger-7.2.0
FEATURE_REF=/wynton/home/pollenlab/jding/BrainChromatin/Perturb/cellranger/feature_ref.csv

SGE_TASK_ID=`expr $SGE_TASK_ID - 1`

SAMPLES=(2D_L2 2D_L3 2D_L4)
#SAMPLES=(HM2D_2nd_L1 HM2D_2nd_L2)
#SAMPLES=(HM2D_L1 HM2D_L2 HM2D_L3)
SAMPLE="${SAMPLES[$SGE_TASK_ID]}"
LIB=/wynton/home/pollenlab/jding/BrainChromatin/Perturb/cellranger/$SAMPLE.csv

echo $SGE_TASK_ID
echo "Running cellranger"
echo "Job directory: $JOB_DIR "
echo "Sample: $SAMPLE "
echo "Reference: $REF "
echo "Feaure barcode: $FEATURE_REF "
echo "Library file: $LIB "



mkdir $JOB_DIR
cd $JOB_DIR


/wynton/group/pollen/jding/cellranger/cellranger-7.2.0/cellranger count --id=$SAMPLE \
           --libraries=$LIB \
           --transcriptome=$REF \
           --feature-ref=$FEATURE_REF \
           #--create-bam true \
           --localcores=4 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"





#! /bin/bash
#$ -e ./log
#$ -o ./log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=50G
#$ -l h_rt=1:00:00
#$ -l scratch=100G
#$ -t 1-3:1
           

source ~/.bashrc
conda activate cellbouncer

JOB_DIR=/wynton/scratch/jding/cellbouncer

SGE_TASK_ID=`expr $SGE_TASK_ID - 1`

SAMPLES=(2D_L2 2D_L3 2D_L4)
#SAMPLES=(HM2D_2nd_L1 HM2D_2nd_L2)
#SAMPLES=(HM2D_L1 HM2D_L2 HM2D_L3)

SAMPLE="${SAMPLES[$SGE_TASK_ID]}"

H5=/wynton/group/pollen/jding/brainchromatin/perturb/cellranger/$SAMPLE/outs/filtered_feature_bc_matrix.h5


echo $SGE_TASK_ID
echo "Running cellbouncer"
echo "Job directory: $JOB_DIR "
echo "Sample: $SAMPLE "
echo "h5 file: $H5 "

mkdir $JOB_DIR
mkdir $JOB_DIR/$SAMPLE
cd $JOB_DIR/$SAMPLE


#convert cellranger h5 to mex
/wynton/home/pollenlab/jding/cellbouncer/utils/h5tomex.py -H $H5 -o $JOB_DIR/$SAMPLE --feature_type "CRISPR Guide Capture" 

#sgRNA assignment
/wynton/home/pollenlab/jding/cellbouncer/demux_tags --sgRNA --cell_barcodes $JOB_DIR/$SAMPLE/barcodes.tsv.gz --features $JOB_DIR/$SAMPLE/features.tsv.gz --mtx $JOB_DIR/$SAMPLE/matrix.mtx.gz --feature_type "CRISPR Guide Capture" -o $JOB_DIR/$SAMPLE


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"





#! /bin/bash
#$ -e ./log2
#$ -o ./log2
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=64G
#$ -l h_rt=5:00:00
#$ -l scratch=200G
#$ -q gpu.q
#$ -l gpu_mem=30000M
#$ -t 1-5:1

#removing ambient RNA using cellbender

source ~/.bashrc

JOB_DIR=/wynton/scratch/jding/Cellranger-7.2.0

SGE_TASK_ID=`expr $SGE_TASK_ID - 1`

SAMPLES=(2D_L2 2D_L3 2D_L4)
#SAMPLES=(HM2D_2nd_L1 HM2D_2nd_L2)
#SAMPLES=(HM2D_L1 HM2D_L2 HM2D_L3)
SAMPLE="${SAMPLES[$SGE_TASK_ID]}"

#run cellranger-7.2.0
cd $JOB_DIR


#run cellbenderv3.0.0
conda activate Cellbender
module load cuda/10.1

export CUDA_VISIBLE_DEVICES=$SGE_GPU

EPOCH=150
MAXCELLS=50000
ADDCELLS=20000


OUT_DIR=/wynton/scratch/jding/Cellbender
mkdir $OUT_DIR


#expected cell number
BARCODE=$JOB_DIR/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
FILE=$JOB_DIR/$SAMPLE/outs/raw_feature_bc_matrix.h5

echo "Running cellbender"
echo "Job directory: $JOB_DIR "
echo "Sample: $SAMPLE "
echo "Input file: $FILE "
echo "Barcode file: $BARCODE "


mkdir $OUT_DIR/$SAMPLE
cd $OUT_DIR/$SAMPLE

CELLS=$(zcat ${BARCODE} | wc -l)
if [ -z "$CELLS" ];then CELLS=$MAXCELLS ; fi
MINCELLS=$(($CELLS < $MAXCELLS ? $CELLS:$MAXCELLS))
MINCELLS=$(($CELLS > 10 ? $CELLS:$MAXCELLS))
DROPS=`expr $MINCELLS + $ADDCELLS`

echo "Expected-cells: $MINCELLS "

cellbender remove-background \
                 --cuda \
                 --input $FILE \
                 --output $OUT_DIR/$SAMPLE/cellbended.h5 \
                 --expected-cells $MINCELLS \
                 --epochs $EPOCH \
                 --learning-rate 1e-4 


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"








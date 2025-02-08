#!/bin/bash
#$ -e ./log
#$ -o ./log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=64G
#$ -l h_rt=60:00:00
#$ -l scratch=120G
#$ -t 2-3:1

source ~/.bashrc
conda activate velocyto

REPEATS=/wynton/home/pollenlab/jding/hg38_rmsk.gtf
GENES=/wynton/group/pollen/genomes/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
JOB_DIR=/wynton/group/pollen/jding/brainchromatin/HM2D/cellranger-7.2.0
OUT_DIR=/wynton/scratch/jding/velocyto

SGE_TASK_ID=`expr $SGE_TASK_ID - 1`
SAMPLES=($((ls $JOB_DIR/*/outs/raw_feature_bc_matrix.h5) | cut -d'/' -f 9))
SAMPLE="${SAMPLES[$SGE_TASK_ID]}"

#use filtered barcodes
BARCODE=$JOB_DIR/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
BAM=$JOB_DIR/$SAMPLE/outs/possorted_genome_bam.bam

echo $REPEATS
echo $GENES
echo $BARCODE
echo $BAM
echo $SAMPLE
echo $OUT_DIR/$SAMPLE


mkdir $OUT_DIR
mkdir $OUT_DIR/$SAMPLE
cd $OUT_DIR/$SAMPLE

                

echo "#####################Making sorted BAM####################"
module load CBI samtools
samtools index $OUT_DIR/$SAMPLE/possorted_genome_bam.bam
samtools sort -t CB -O BAM -o $OUT_DIR/$SAMPLE/cellsorted_possorted_genome_bam.bam $BAM

#echo "########################Copying GTF###########################"
cp $GENES $OUT_DIR/$SAMPLE
gunzip $OUT_DIR/$SAMPLE/genes.gtf.gz

echo "#####################Running Velocyto#######################"
velocyto run -b $BARCODE -o $OUT_DIR/$SAMPLE/ -m $REPEATS $OUT_DIR/$SAMPLE/possorted_genome_bam.bam \
             $OUT_DIR/$SAMPLE/genes.gtf \
              --samtools-memory 640000
             
#rm $OUT_DIR/$SAMPLE/cellsorted_possorted_genome_bam.bam
#rm $OUT_DIR/$SAMPLE/possorted_genome_bam.bam
rm $OUT_DIR/$SAMPLE/genes.gtf

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


#! /bin/bash
#$ -e ./log3
#$ -o ./log3
#$ -cwd
#$ -j y
#$ -pe mpi 16
#$ -l mem_free=100G
#$ -l h_rt=80:00:00
#$ -l scratch=200G
#$ -t 1-3:1

#running vireo to assigning donors to cells based on SNP genotypes

source ~/.bashrc
conda activate vireo

SGE_TASK_ID=`expr $SGE_TASK_ID - 1`


REGION_VCF=/wynton/group/pollen/jding/vireo/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz
OUT_DIR=/wynton/scratch/jding/Vireo

JOB_DIR=/wynton/scratch/jding/Cellranger-7.2.0
FILES=($JOB_DIR/*/outs/possorted_genome_bam.bam)
BAM="${FILES[$SGE_TASK_ID]}"
BARCODES=($JOB_DIR/*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz)
BARCODE="${BARCODES[$SGE_TASK_ID]}"

SAMPLES=($((ls $JOB_DIR/*/outs/possorted_genome_bam.bam) | cut -d'/' -f 6))
SAMPLE="${SAMPLES[$SGE_TASK_ID]}"

OUT=$OUT_DIR/$SAMPLE

n_donor=3

echo $JOB_DIR
echo $SAMPLE
echo $BARCODE
echo $BAM
echo $REGION_VCF
echo $OUT


mkdir $OUT_DIR
cd $OUT_DIR
mkdir $OUT

cellsnp-lite -s $BAM -b $BARCODE -O $OUT -R $REGION_VCF -p 16 --minMAF 0.1 --minCOUNT 20 --gzip --genotype
vireo -c $OUT -N $n_donor -o $OUT -p 16

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"



           


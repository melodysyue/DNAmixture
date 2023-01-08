#!/bin/bash -l

#SBATCH --job-name=bwa_ntc
#SBATCH -o bwa_ntc.out
#SBATCH -e bwa_ntc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/DNAmixture/data/ntc/trimSE


enable_lmod
module load bwa
module load samtools

files=$(find . -name '*_L001_trim.fastq.gz' | sort )

for file in $files
do 
   OUT="$(basename $file _L001_trim.fastq.gz)"
   bwa mem -t 31 ../../../reference/CHreference_oneline.fasta $file | \
   samtools sort -@31 -o ${OUT}_sorted.bam -;
done;


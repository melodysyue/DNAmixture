#!/bin/bash -l

#SBATCH --job-name=samIndex_ntc
#SBATCH -o samIndex_ntc.out
#SBATCH -e samIndex_ntc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/DNAmixture/data/ntc/trimSE


enable_lmod
module load samtools

bam=$(find . -name '*_sorted.bam' | sort )

for file in $bam
do 
   samtools index $file;
done;


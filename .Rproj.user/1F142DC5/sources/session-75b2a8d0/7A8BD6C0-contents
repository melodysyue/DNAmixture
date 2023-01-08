#!/bin/bash -l

#SBATCH --job-name=mapstats_ntc
#SBATCH -o mapstats_ntc.out
#SBATCH -e mapstats_ntc.err
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
   idx="$(basename $file _sorted.bam)_idxstats.txt"
   flag="$(basename $file _sorted.bam)_flagstats.txt"
   samtools idxstats $file > $idx
   samtools flagstat $file > $flag;
done


#!/bin/bash -l

#SBATCH --job-name=b2s_ccf
#SBATCH -o b2s_ccf.out
#SBATCH -e b2s_ccf.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/DNAmixture/data/ccf_fin/trimSE


enable_lmod
module load samtools

bam=$(find . -name '*_sorted.bam' | sort )

for file in $bam
do 
   sam="$(basename $file _sorted.bam).sam"
   samtools view -h $file > $sam;
done


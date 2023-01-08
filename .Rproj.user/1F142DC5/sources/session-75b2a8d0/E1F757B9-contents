#!/bin/bash -l

#SBATCH --job-name=trimSE_ntc
#SBATCH -o trimSE_ntc.out
#SBATCH -e trimSE_ntc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 31
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish
#SBATCH --chdir=/scratch-lustre/wlarson/yue/DNAmixture/data/ntc


enable_lmod
module load trimmomatic
#module load java

mkdir trimSE

forward=$(find . -name '*_R1_001.fastq.gz' | sort )

for R1 in $forward
do 
   LOG="$(basename $R1 _R1_001.fastq.gz).log"
   OUT="$(basename $R1 _R1_001.fastq.gz)_trim.fastq.gz"
   trimmomatic SE -trimlog ./trimSE/$LOG -threads 31 $R1 ./trimSE/$OUT ILLUMINACLIP:../../script/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50 
done

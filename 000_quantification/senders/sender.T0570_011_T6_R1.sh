#!/bin/bash

#SBATCH --job-name=T0570_011_T6_R1
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.T0570_011_T6_R1.out.txt   
#SBATCH --error=messages/messages.T0570_011_T6_R1.err.txt 

date
pwd

#
# 1. Create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users
if [ ! -d $scratchlocation/$USER ]; then
mkdir -p $scratchlocation/$USER
fi
tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

#
# 2. copy files into scratch
#
echo ""
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/akthelia/data/T0570_011_T6_R1 $tdir/.
date


#
# 3. call Trimmomatic
#
echo ""
echo "about to call trimmomatic"
date
cd $tdir
mkdir $tdir/clean_fastq
time java -jar /users/home/adrian/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 T0570_011_T6_R1/T0570_011_T6_R1_1.fq.gz T0570_011_T6_R1/T0570_011_T6_R1_2.fq.gz clean_fastq/T0570_011_T6__R1_clean.fastq.gz clean_fastq/T0570_011_T6__R1_garbage.fastq.gz clean_fastq/T0570_011_T6__R2_clean.fastq.gz clean_fastq/T0570_011_T6__R2_garbage.fastq.gz ILLUMINACLIP:/users/home/adrian/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_a -t 4 -b 100 --rf-stranded --verbose clean_fastq/T0570_011_T6__R1_clean.fastq.gz clean_fastq/T0570_011_T6__R2_clean.fastq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_b -t 4 -b 100 --fr-stranded --verbose clean_fastq/T0570_011_T6__R1_clean.fastq.gz clean_fastq/T0570_011_T6__R2_clean.fastq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_c -t 4 -b 100 --verbose clean_fastq/T0570_011_T6__R1_clean.fastq.gz clean_fastq/T0570_011_T6__R2_clean.fastq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
date
mkdir /hpcdata/Mimir/adrian/research/akthelia/results/T0570_011_T6_R1_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/akthelia/results/T0570_011_T6_R1_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    
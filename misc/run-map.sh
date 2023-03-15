#!/bin/bash
#SBATCH -J aln
#SBATCH -p uohhm
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=7000
#SBATCH --mail-user=digenova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o aln_%j.out
#SBATCH -e aln_%j.err
#sleep 60
#echo "Hosthame: `hostname`"


#iodd.py --metaplasmid  -1 ./reads/17.R1.fq.gz -2 ./reads/17.R2.fq.gz -t 8 -m 80 -o 17_iodd_${SLURM_ARRAY_TASK_ID}

make all 

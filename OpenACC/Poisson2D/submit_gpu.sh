#! /bin/bash
#
## Walltime limit
#$ -l h_rt=0:05:00
#
## Give the job a name.
#$ -N acc_test
#
## Redirect error output to standard output
#$ -j y
#
## What project to use. "paralg" is the project for the class
#$ -P ec526
#
## Ask for nodes with 4 cores
#$ -pe omp 4
#
## Ask for one GPU with at least compute capability 5.0
#$ -l gpus=0.25
#$ -l gpu_c=5.0

# Want more flags? Look here:
# http://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/


# Immediately form fused output/error file, besides the one with the default name.
exec >  ${SGE_O_WORKDIR}/${JOB_NAME}-${JOB_ID}.scc.out 2>&1

./poisson2d

exit


#!/bin/bash -l
# Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cores per node required 
#SBATCH -n 44

# Set the walltime of the job to 1 hour (format is hh:mm:ss)
#SBATCH -t 168:00:00

# E-mail on begin (b), abort (a) and end (e) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=barry.wardell@ucd.ie

# Specifies the jobname
#SBATCH --job-name=d2R-SS-r7

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR

echo "Running on `hostname`"
echo "Started at `date`"

module load mathematica
WolframKernel -noprompt -run "r0=7.0;" -script d2R-SS-r0.wl

echo "Done at `date`"

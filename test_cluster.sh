#!/bin/bash

#SBATCH --job-name=testjob
#SBATCH --chdir=/home/

# Resources
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00-08:00:00

# User config
#SBATCH --mail-type=FAIL,SUCCESS
#SBATCH --mail-user=amanda.nicotina@tum.de
#SBATCH --output /bnmrz/sg/nicotina/Slurm/Documents/slurm-%j.out
#SBATCH --error /bnmrz/sg/nicotina/Slurm/Documents/slurm-%j.err

# Runtime

cd /bnmrz/sg/nicotina/Documents/Slurm/project1/ || exit 1

julia test_script.jl || exit 2
#!/usr/bin/env bash
#SBATCH --job-name=
#SBATCH --time 5-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH -o %j.out
#SBATCH -e %j.err


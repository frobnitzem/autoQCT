#!/bin/bash
#SBATCH --nodes=3
#SBATCH --time=120
#SBATCH --job-name=autoQCT
#SBATCH --output=%j.out

export PYTHONPATH=$PYTHONPATH:$HOME/autoQCT
export PSI_SCRATCH=/qscratch/$USER

PMAKE=$HOME/pmake/make.py
$PMAKE rules.yaml targets.yaml 118


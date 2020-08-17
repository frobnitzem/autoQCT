#!/bin/bash
#SBATCH --nodes=10
#SBATCH --time=45:00
#SBATCH --job-name=autoQCT
#SBATCH --output=%j.out

export PYTHONPATH=$PYTHONPATH:$HOME/autoQCT
PMAKE=$HOME/pmake/make.py
$PMAKE rules.yaml targets.yaml 45


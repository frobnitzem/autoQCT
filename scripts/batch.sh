#!/bin/bash
#SBATCH --nodes=10
#SBATCH --time=10:00
#SBATCH --job-name=autoQCT
#SBATCH --output=%j.out

python pmake.py rules.yaml targets.yaml


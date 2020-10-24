# -*- coding: utf-8 -*-

""" Creates slurm jobs batching a list of QM calculations. """

__all__ = [ 'mk_slurm' ]

scr = """#!/bin/bash
#SBATCH -J BatchGaussian
%(sbt)s

module load gaussian09

mpirun --bynode -np 1 %(procs)s
"""


def mk_slurm(coms, **args):
    # Create the text of a slurm job file, using keyword-arguments
    # to set the the "#SBATCH --key=value" list.
    # 
    # Note: This replaces all underscores with dashes in the keys.
    #
    # Inputs:
    #
    #   coms = list of strings, each the name of a gaussian com-file
    #   **args = slurm parameters
    #
    # Output:
    #
    #   Returns a string with the slurm job contents.
    #
    # Example Usage:
    #   mk_slurm(["1.com", "2.com", "3.com"],
    #       partition='defq',
    #       qos='normal',
    #       time='1:00:00',
    #       nodes='1',
    #       ntasks_per_node='1',
    #       cpus_per_task='20',
    #       gres='mic:0',
    #       mem='64000',
    #       export='ALL')
    #

    sbt = []
    for k, v in args.items():
        k = k.replace("_", "-")
        sbt.append("#SBATCH --%s=%s"%(k, v))

    return scr % {'sbt': '\n'.join(sbt),
                  'procs': ' : '.join("g09 %s"%x for x in coms)
                 }


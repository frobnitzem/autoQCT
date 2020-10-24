# -*- coding: utf-8 -*-

""" Creates QM input files for running energy calculations. """

import numpy

__all__ = [ 'mk_com' ]

scr = """%%chk=%(name)s.chk
%%NProcShared=%(cores)d
%(route)s

%(name)s

%(charge)d %(spin)d
%(crds)s
"""

def mk_com(mol, name, cores, route):
    # Create the text of a Gaussian09 / 16 input file.
    #
    # Inputs:
    #
    #   mol   : Mol/Sys = molecule or system to calculate on
    #   name  : string = unique file name
    #   cores : int = number of CPU cores to use
    #   route : string = Gaussian "route"
    #
    # Output:
    #
    #   Returns a string with the 'com'-file contents.
    #
    # Example Usage:
    #   wat = Mol(["O", "H", "H"], [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]], 0, 1)
    #   mk_com(wat, "wat_en", 16, "#N B3LYP/aug-cc-pvdz SP NOSYMMETRY")

    crds = []
    for name, x in zip(mol.name_iter(), mol.x_iter()):
        crds.append("%4s  %12.6f %12.6f %12.6f"%(name,x[0],x[1],x[2]))

    return scr % {
            'name': name,
            'cores': cores,
            'route': route,
            'charge': mol.tot_charge(),
            'spin': mol.tot_spin(),
            'crds': '\n'.join(crds)
           }


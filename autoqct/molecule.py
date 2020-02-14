# -*- coding: utf-8 -*-
  
"""Definitions for chemicals."""

__all__ = [ 'Mol', 'Sys' ]

import numpy

class Mol:
    def __init__(self, names, x, charge=0, spin=1):
        # Inputs:
        #   names : [string] = element names for all N atoms
        #   x     : array, shape=(N,3), dtype=float = atomic coordinates
        #   charge : int = charge
        #   spin   : int = (total electron spin)*2 + 1
        #
        # Note:
        #
        #   The total electron spin is the number of unpaired
        #   electrons times 1/2.
        #
        self.names = names
        self.x = x
        self.charge = charge
        self.spin = spin

        assert isinstance(x, numpy.ndarray)
        assert len(names) == x.shape[0]
        assert x.shape[1] == 3
        assert len(x.shape) == 2

    def name_iter(self): # return an iterator through all atom names
        return self.names

    def x_iter(self): # return an iterator through all atom crds
        return self.x

    def tot_charge(self):
        return self.charge
    def tot_spin(self):
        return self.spin

class Sys:
    """ A list of molecules. """

    def __init__(self, mols):
        self.mols = mols
    def name_iter(self):
        for m in self.mols:
            for a in m.name_iter():
                yield a
    def x_iter(self):
        for m in self.mols:
            for x in m.x_iter():
                yield x

    def tot_charge(self):
        return sum(m.tot_charge() for m in self.mols)
    def tot_spin(self):
        return sum(m.tot_spin()-1 for m in self.mols)+1


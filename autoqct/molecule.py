# -*- coding: utf-8 -*-
  
"""Definitions for chemicals."""

__all__ = [ 'Mol', 'Sys' ]

class Mol:
    def __init__(self, names, charge=0, spin=1):
        # Inputs:
        #   names : [string] = element names for all N atoms
        #   charge : int = charge
        #   spin   : int = (total electron spin)*2 + 1
        #
        # Note:
        #
        #   The total electron spin is the number of unpaired
        #   electrons times 1/2.
        #
        self.names = names
        self.__atoms = len(names)
        self.__charge = charge
        self.__spin = spin

    def copy(self):
        return Mol(self.names, self.charge, self.spin)

    def name_iter(self): # return an iterator through all atom names
        return self.names

    def atoms(self):
        return self.__atoms
    def charge(self):
        return self.__charge
    def spin(self):
        return self.__spin

class Sys(list):
    """ A list of molecules. """

    def __init__(self, mols=[]):
        list.__init__(self, mols)

    def name_iter(self):
        for m in self:
            for a in m.name_iter():
                yield a
    def atoms(self):
        return sum(m.atoms() for m in self)
    def charge(self):
        return sum(m.charge() for m in self)
    def spin(self):
        return sum(m.spin()-1 for m in self)+1

    def fmt_psi4(self, x, detect_symm=True):
        assert len(x) == self.atoms()
        s = [ "%d %d"%(self.charge(), self.spin()) ]
        if detect_symm == False:
            s.append("symmetry c1")
        fmt = "%2s %.5f %.5f %.5f"
        s.extend( fmt%(n,y[0],y[1],y[2]) for n,y in zip(self.name_iter(), x) )
        return "\n".join(s)


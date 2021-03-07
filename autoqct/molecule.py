# -*- coding: utf-8 -*-
  
"""Definitions for chemicals."""

__all__ = [ 'Mol', 'Sys','elems','elem_dict' ]

elems = [] # (name, mass)
elem_dict = {}
for i, line in enumerate("""0 Unk 1.0
1   H 1.008
2   He 4.0026
3   Li 6.94
4   Be 9.0122
5   B 10.81
6   C 12.011
7   N 14.007
8   O 15.999
9   F 18.998
10  Ne 20.180
11  Na 22.990
12  Mg 24.305
13  Al 26.982
14  Si 28.085
15  P 30.974
16  S 32.06
17  Cl 35.45
18  Ar 39.948
19  K 39.098
20  Ca 40.078
21  Sc 44.956
22  Ti 47.867
23  V 50.942
24  Cr 51.996
25  Mn 54.938
26  Fe 55.845
27  Co 58.933
28  Ni 58.693
29  Cu 63.546
30  Zn 65.38""".split('\n')):
    t = line.split()
    assert int(t[0]) == i
    elems.append((t[1], float(t[2])))
    elem_dict[t[1]] = (i, float(t[2]))

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


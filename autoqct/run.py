# This script intended to be run by extending the QCT
# class and then calling autoQCT(Q)
#
# autoQCT reads sys.argv in the following formats:
#
# python {inp[qct]} ebind {inp[crd]} {out[en]}
# python {inp[qct]} emin data/frame.12.npy {out[en]}
# python {inp[qct]} exp_pos nm1_%d.en n_%d.en {frames} {out[dg]}
# python {inp[qct]} solv data/frame.12.npy data/solv.12.json
# python {inp[qct]} exp_neg n_%d.muex {frames} {out[solv]}
# python {inp[qct]} rrho data/frame.12.npy rrho.yaml [0-based molecule numbers to minimize]

import os, json
import numpy as np
from pathlib import Path

from .molecule import *
from .cp2kfile import read_cp2k
from .pcm import run_pcm
from .ebind import run_ebind
from .rrho import run_rrho, run_emin
from .importer import import_py

__all__ = ['QCT', 'autoQCT', 'dist']

kT  = 8.3145e-3*298.15 / 4.184 # kcal/mol

# helper function - calculate pair distance
def dist(x, y):
    return np.sqrt( np.sum((x-y)*(x-y)) )

# helper function - read all json files into a list
def read_json(name_fmt, frames, key=None):
    en = []
    for i in range(frames):
        with open(name_fmt % i) as f:
            e = json.loads(f.read())
            if key is not None:
                e = e[key]
            en.append( e )
    return en

class QCT:
    attr = {
        'molecules': list,
        'clustered': 'function',
        'theory':    str,
        'basis':     str
    }
    def __init__(self):
        self.name = self.__class__.__name__
        for k, v in self.attr.items():
            if not hasattr(self, k):
                raise NameError("Incomplete user code: '%s' is not defined."%k)
        self.system = Sys([m[1] for m in self.molecules])

    # Test the filter by showing the frame 'fr' and its run-result.
    def test_filter(self, trj_name, fr="0"):
        trj = read_cp2k(trj_name, self.molecules)

        i = int(fr)
        if len(trj) <= i:
            print("Missing frame")
            return
        x = trj[i].x
        print( self.clustered(x) )
        r1 = np.sum( (x[1::3]-x[2::3])**2, -1 )
        r2 = np.sum( (x[1::3]-x[3::3])**2, -1 )
        if np.any(r1 > 1.1) or np.any(r2 > 1.1):
            print("bond distances out of bounds:")
            print(np.sqrt(r1))
            print(np.sqrt(r2))
            print(x)
        else:
            print("bond distances OK")

    # TODO: turn read_cp2k into an iterator
    def init(self, trj_name, start="-1"):
        trj = read_cp2k(trj_name, self.molecules)
        start = int(start)
        if start == -1:
            start = len(trj)//2
        trj = trj[start:]
        #for i in range(start-1, -1, -1):
        #    del trj[i]

        print("Read trajectory with %d frames." % len(trj))
        for i in range(len(trj)-1, -1, -1):
            if not self.clustered(trj[i].x):
                del trj[i]

        print("%d frames pass filter" % len(trj))

        targets = """{name}:
    frames: {frames:d}
    start:  {start:d}
    dirname: data
    out:
        dg: dg_nm1_n.en
        solv: dg_solv.json
        thermo: nm1.json
""".format(name=self.name, frames=len(trj), start=start)
        with open("targets.yaml", "w") as f:
            f.write(targets)
        data = Path("data")
        data.mkdir(exist_ok=False)
        #data.mkdir(exist_ok=True)
        for i,frame in enumerate(trj):
            np.save(data / ("frame_%d.npy"%i), frame.x)

    def reduce(self, trj_name, fn_path):
        module = import_py("user_fn", fn_path)
        trj = read_cp2k(trj_name, self.molecules)

        #for i in range(len(trj)-1, -1, -1):
        #    if not self.clustered(trj[i].x):
        #        del trj[i]
        x = module.x0
        for i in range(len(trj)):
            module.fn(self, trj[i].x, x)
        module.output(trj_name, x)

    def loop_X_nm1(self, x):
        # Loop over n solvent molecules in Xn and make pairs
        # of X1 and X(n-1) systems
        n0 = self.system[0].atoms()
        k = n0
        de = np.zeros(len(self.system)-1)
        N = self.system.atoms()
        for j in range(1, len(self.system)):
            n = self.system[j].atoms()

            # splice out mol. j
            s1 = Sys( self.system[:j] + self.system[j+1:] )
            x1 = np.empty((N-n,3))
            x1[:k] = x[:k]
            x1[k:] = x[k+n:]

            # splice X with mol j
            s2 = Sys( self.system[:1] + self.system[j:j+1] )
            x2 = np.empty((n0+n,3))
            x2[:n0] = x[:n0]
            x2[n0:] = x[k:k+n]

            yield s1, x1, s2, x2
            k += n

    def loop_nm1(self, x):
        # Loop over solvent molecules and make "leave-one-out"
        # nm1 = n - 1 size clusters.
        k = self.system[0].atoms()
        de = np.zeros(len(self.system)-1)
        N = self.system.atoms()
        for j in range(1, len(self.system)):
            n = self.system[j].atoms()

            # splice out mol. j
            s1 = Sys( self.system[:j] + self.system[j+1:] )
            x1 = np.empty((N-n,3))
            x1[:k] = x[:k]
            x1[k:] = x[k+n:]

            yield s1, x1, Sys( [self.system[j]] ), x[k:k+n]

            k += n

    def ecluster(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        k = self.system[0].atoms() # number of atoms in solute
        x1 = x[:k]
        x2 = x[k:]
        s1 = Sys( self.system[:1] ) # molecule (solute)
        s2 = Sys( self.system[1:] ) # all other molecules
        de = run_ebind(s1, x1, s2, x2, theory=self.theory, basis=self.basis)

        with open(out, "w") as f:
            f.write( "%g\n"%de )

    def ebind_old(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        de = []
        for s1, x1, s2, x2 in self.loop_nm1(x):
            #print(x1)
            #print(x2)
            de.append( run_ebind(s1, x1, s2, x2,
                                 theory=self.theory, basis=self.basis)
                     )
            #print(de[-1])
            #print()

        with open(out, "w") as f:
            f.write( json.dumps(de) + '\n' )

    def ebind(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        e0 = run_emin(self.system, x,
                      theory=self.theory, basis=self.basis,
                      optimize=False)
        de = {'Xn': e0, 'Xnm1': [], 'X1': []}
        for s1, x1, s2, x2 in self.loop_X_nm1(x):
            e1 = run_emin(s1, x1,
                          theory=self.theory, basis=self.basis,
                          optimize=False)
            e2 = run_emin(s2, x2,
                          theory=self.theory, basis=self.basis,
                          optimize=False)
            de['Xnm1'].append(e1)
            de['X1'].append(e2)

        with open(out, "w") as f:
            f.write( json.dumps(de) + '\n' )

    def solv(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)
        ret = run_pcm(self.system, x, theory=self.theory, basis=self.basis)
        with open(out, "w") as f:
            f.write( json.dumps(ret, indent=4) + '\n' )

    # run energy minimization on the system
    def emin(self, xyz, out):
        x = np.load(xyz)
        en = run_emin(self.system, x,
                      theory=self.theory, basis=self.basis)
        if en is None:
            return 1
        with open(out, "w") as f:
            f.write("%f\n"%en)
        return 0

    # run rrho on K lowest energy structures
    def rrho1(self, trj_name, out, K=1, start="-1"):
        trj = read_cp2k(trj_name, self.molecules)
        K = int(K)
        start = int(start)
        if start == -1:
            start = len(trj)//2

        trj = trj[start:]
        E = [(f.energy,i) for i,f in enumerate(trj)]
        E.sort()
        ans = []
        for en,i in E:
            if len(ans) >= K:
                break
            r = run_rrho(self.system, trj[i].x,
                         theory=self.theory, basis=self.basis)
            if r is None:
                continue
            r['init'] = en
            ans.append(r)
            with open(out, "w") as f:
                f.write( json.dumps(ans, indent=4) + '\n' )

    # run rrho on all mols individually
    def rrho(self, xyz, out, *mols):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        ans = []
        for si in mols:
            i = int(si)
            start = sum(s.atoms() for s in self.system[:i])
            k = self.system[i].atoms()
            assert i >= 0 and i < len(self.system)
            ans.append(
                run_rrho(Sys( [self.system[i]] ),
                         x[start:start+k]
                         theory=self.theory, basis=self.basis)
                )
        with open(out, "w") as f:
            f.write( json.dumps(ans, indent=4) + '\n' )

    # en_nm1 and en_n are input filenames
    # dg is an output filename
    def exp_pos(self, name_fmt, nfr, dg):
        frames = int(nfr)
        en = read_json(name_fmt, frames)
        en  = np.array(en).reshape(-1) / kT
        top = en.max()
        avg = np.sum( np.exp(en - top) / len(en) )
        en = kT * (top + np.log(avg))
        with open(dg, "w") as f:
            f.write("%f\n"%en)

    # muex_n is the input filename pattern
    # nfr is the number of frames (0, 1, 2, ..., nfr-1)
    # out is the output filename
    # example: exp_neg solv_%d.json {frames} {out[solv]}
    def exp_neg(self, muex_n, nfr, out):
        frames = int(nfr)
        en = read_json(muex_n, frames, "PCM Polarization")
        en = np.array(en).reshape(-1)
        emin = en.min()
        avg = np.sum(np.exp( -(en - emin)/kT )) / len(en)
        en = emin - kT*np.log( avg )
        with open(out, "w") as f:
            f.write("%f\n"%en)

def autoQCT(Q):
    import sys
    argv = sys.argv
    assert len(argv) >= 3, "Usage: %s <method name> ..."

    q = Q()

    #import importlib
    #mod = importlib.import_module(argv[1])
    # q = mod.system()
    try:
        fn = getattr(q, argv[1])
    except AttributeError:
        print("Error: method '%s' does not exist." % argv[1])
        return 1

    return fn(*argv[2:])


# This script intended to be run by extending the QCT
# class and then calling autoQCT(Q)
#
# autoQCT reads sys.argv in the following formats:
#
# python {inp[qct]} ebind {inp[crd]} {out[en]}
# python {inp[qct]} exp_pos nm1_%d.en n_%d.en {frames} {out[dg]} {out[min]}
# python {inp[qct]} solv {inp[crd]} {out[en]}
# python {inp[qct]} exp_neg n_%d.muex {frames} {out[solv]}
# python {inp[qct]} rrho {inp[min]} {out[thermo]} 

import os, json
import numpy as np
from pathlib import Path

from .molecule import *
from .cp2kfile import read_cp2k
from .pcm import run_pcm
from .ebind import run_ebind
from .rrho import run_rrho

__all__ = ['QCT', 'autoQCT', 'dist']

kT  = 8.3145e-3*298.15 / 4.184

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
        print(trj[i].x)
        print( self.clustered(trj[i].x) )

    # TODO: turn read_cp2k into an iterator
    def init(self, trj_name):
        trj = read_cp2k(trj_name, self.molecules)

        print("Read trajectory with %d frames." % len(trj))
        for i in range(len(trj)-1, -1, -1):
            if not self.clustered(trj[i].x):
                del trj[i]

        print("%d frames pass filter" % len(trj))

        targets = """{name}:
    frames: {frames:d}
    dirname: data
    out:
        dg: dg_nm1_n.en
        solv: dg_solv.json
        thermo: nm1.json
""".format(name=self.name, frames=len(trj))
        with open("targets.yaml", "w") as f:
            f.write(targets)
        data = Path("data")
        data.mkdir(exist_ok=True)
        for i,frame in enumerate(trj):
            np.save(data / ("frame_%d.npy"%i), frame.x)

    def loop_nm1(self, x):
        # Loop over solvent molecules and make "leave-out-out"
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
        s1 = Sys( self.system[:1] )
        s2 = Sys( self.system[1:] )
        de = run_ebind(s1, x1, s2, x2, theory=self.theory, basis=self.basis)

        with open(out, "w") as f:
            f.write( "%g\n"%de )

    def ebind(self, xyz, out):
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

    def solv(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)
        ret = run_pcm(self.system, x, theory=self.theory, basis=self.basis)
        with open(out, "w") as f:
            f.write( json.dumps(ret, indent=4) + '\n' )

    # run rrho on all nm1 clusters
    def rrho(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        ans = {}
        ans['solute'] = run_rrho(Sys( [ self.system[0] ] ),
                                 x[:self.system[0].atoms()],
                                 theory=self.theory, basis=self.basis)
        ans['solvent'] = {}
        ans['nm1'] = []

        for s1, x1, s2, x2 in self.loop_nm1(x):
            ret = run_rrho(s1, x1, theory=self.theory, basis=self.basis)
            ans['nm1'].append(ret)
            # FIXME: use better solvent naming scheme
            solv_name = ''.join(s2.name_iter()) # e.g. 'OHH'
            if solv_name not in ans['solvent']:
                ans['solvent'][solv_name] = \
                    run_rrho(s2, x2, theory=self.theory, basis=self.basis)
        with open(out, "w") as f:
            f.write( json.dumps(ans, indent=4) + '\n' )

    # en_nm1 and en_n are input filenames
    # dg, minxyz are output filenames
    def exp_pos(self, name_fmt, nfr, dg, minxyz):
        frames = int(nfr)
        en = read_json(name_fmt, frames)
        en  = np.array(en).reshape(-1) / kT
        m = np.argmax(en) # use worst binder -- maybe n-1 is "already good"
        top = en.max()
        avg = np.sum( np.exp(en - top) / len(en) )
        en = kT * (top + np.log(avg))
        with open(dg, "w") as f:
            f.write("%f\n"%en)
        import shutil
        shutil.copy("frame_%d.npy"%m, minxyz)

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


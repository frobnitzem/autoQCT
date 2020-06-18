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
from .cp2kfile import *
from .qm_input import *

__all__ = ['QCT', 'autoQCT']

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

    def init(self, trj_name):
        trj = read_cp2k(trj_name, self.molecules)

        print("Read trajectory with %d frames." % len(trj))
        for i in range(len(trj)-1, -1, -1):
            if not self.clustered(trj[i].x):
                del trj[i]

        print("Filtered out %d frames." % len(trj))

        targets = """{name}:
    frames: {frames:d}
    out:
        dg: dg_nm1_n.en
        solv: dg_solv.txt
        thermo: nm1.yaml
""".format(name=self.name, frames=self.frames)
        with open("targets.yaml", "w") as f:
            f.write(targets)
        for i,frame in enumerate(trj):
            np.save("frame_%d.npy", frame.x)

    def loop_nm1(self, x):
        # Loop over solvent molecules and make "leave-out-out"
        # nm1 = n - 1 size clusters.
        k = self.system[0].atoms()
        de = np.zeros(len(self.system)-1)
        for j in range(1, len(self.system)):
            n = self.system[j].atoms()

            # splice out mol. j
            s1 = Sys( self.system[:j] + self.system[j+1:] )
            x1 = np.empty((N-n,3))
            x1[:k] = x[:k]
            x1[k:] = x[k+n:]

            yield s1, x1, Sys( [self.system[j]] ), x[k:k+n]

            k += n

    def ebind(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        de = []
        for s1, x1, s2, x2 in self.loop_nm1(x):
            de.append( run_ebind(s1, x1, s2, x2,
                                 theory=self.theory, basis=self.basis)
                     )

        with open(out, "w") as f:
            out.write( json.dumps(de) + '\n' )

    def solv(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)
        ret = run_pcm(self.system, x, theory=self.theory, basis=self.basis)
        with open(out, "w") as f:
            out.write( json.dumps(ret, indent=4) + '\n' )

    # run rrho on all nm1 clusters
    def rrho(self, xyz, out):
        x = np.load(xyz)
        N = self.system.atoms()
        assert x.shape == (N,3)

        ans = {}
        ans['solute'] = run_rrho(Sys( self.system[0] ),
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
            out.write( json.dumps(ans, indent=4) + '\n' )

    # en_nm1 and en_n are input filenames
    # dg, minxyz are output filenames
    def exp_pos(self, name_fmt, nfr, dg, minxyz):
        frames = int(nfr)
        en = []
        for i in range(frames):
            with open(name_fmt % i) as f:
                en.append( json.loads(f.read()) )
        kT  = 8.3145e-3*298.15 / 4.184
        en  = np.array(en).reshape(-1) / kT
        top = en.max()
        avg = np.sum( np.exp(en - top) / len(en) )
        return kT * (top + np.log(avg))

    # muex_n is the input filename
    # out is the output filename
    def exp_neg(self, muex_n, nfr, out):
        frames = int(nfr)
        with open(out, "w") as f:
            out.write("0.0\n")

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


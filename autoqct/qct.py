#!/usr/bin/env python3

#python {inp[qct]} singlept {inp[crd]} {out[en]}
#python {inp[qct]} exp_pos nm1_%d.en n_%d.en {frames} {out[dg]} {out[min]}
#python {inp[qct]} solv {inp[crd]} {out[en]}
#python {inp[qct]} exp_neg n_%d.muex {frames} {out[solv]}
#python {inp[qct]} rrho {inp[min]} {out[thermo]} 

import importlib

def singlept(sys, xyz, out):
    with open(out, "w") as f:
        out.write("0.0\n")

def solv(sys, xyz, out):
    with open(out, "w") as f:
        out.write("0.0\n")

def rrho(sys, xyz, out):
    #yaml.write(out)
    pass

# en_nm1 and en_n are input filenames
# dg, minxyz are output filenames
def exp_pos(sys, en_nm1, en_n, nfr, dg, minxyz):
    frames = int(nfr)

# muex_n is the input filename
# out is the output filename
def exp_neg(sys, muex_n, nfr, out):
    frames = int(nfr)
    with open(out, "w") as f:
        out.write("0.0\n")

impl = {'singlept' : singlept,
        'solv' : solv,
        'rrho' : rrho,
        'exp_pos' : exp_pos,
        'exp_neg' : exp_neg
       }

def main(argv):
    assert len(argv) >= 3, "Usage: %s <sys.py> <method name> ..."

    mod = importlib.import_module(argv[1])
    #mod.HelloWorld()
    if argv[2] not in impl:
        print("Error: method %s does not exist." % argv[2])
        return 1

    impl[argv[2]](mod, *argv[3:])
    return 0

if __name__ == "__main__":
    import sys
    exit( main(sys.argv) )


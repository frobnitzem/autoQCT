from tempfile import TemporaryDirectory
from pathlib import Path
from subprocess import check_output, CalledProcessError, STDOUT
import re, os

class Psi4Run:
    def __init__(self, inp, base='/dev/shm'):
        self.inp = inp
        #base = "/qscratch/%s" % os.environ['USER']
        self.tmpdir = TemporaryDirectory(dir=base)
        self.dirname = Path( self.tmpdir.name )
        self.infile = self.dirname / "input.dat"
        self.err = None

    def run(self):
        with open(self.infile, "w") as f:
            f.write(self.inp)
        try:
            #self.out = check_output(["cat input.dat"],
            #             stderr=STDOUT, shell=True,
            #             cwd=str(self.dirname))
            self.out = check_output(["psi4"], #, "-n", "18"],
                         stderr=STDOUT,
                         cwd=str(self.dirname))
            #print(self.out)
        except CalledProcessError as e:
            print("Received error code %d"%e.returncode)
            print(e.output.decode('utf-8'))
            self.err = e.returncode
            self.out = e.output

    # if True, scale all floats by psi4.constants.hartree2kcalmol
    def scrape(self, regex, scale=False):
        hartree2kcalmol = 627.503
        attr = {}
        expr = re.compile(regex)
        fexpr = re.compile(r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?')
        with open(self.dirname / "output.dat") as f:
            for line in f:
                m = expr.match(line)
                if m is None:
                    continue
                v = fexpr.match(m[2])
                if v is None:
                    continue
                s = float(v[0])
                if scale:
                    s *= hartree2kcalmol
                attr[m[1].strip()] = s
        return attr


from tempfile import TemporaryDirectory
from pathlib import Path
from subprocess import check_output, CalledProcessError, STDOUT
import re

class Psi4Run:
    def __init__(self, inp):
        self.inp = inp
        self.tmpdir = TemporaryDirectory()
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
            self.out = check_output(["psi4"],
                         stderr=STDOUT,
                         cwd=str(self.dirname))
            #print(self.out)
        except CalledProcessError as e:
            print("Received error code %d"%e.returncode)
            print(e.output)
            self.err = e.returncode
            self.out = e.output

    def scrape(self, regex):
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
                attr[m[1].strip()] = float(v[0])
        return attr


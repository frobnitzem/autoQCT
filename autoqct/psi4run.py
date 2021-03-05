from tempfile import TemporaryDirectory, mkdtemp
from pathlib import Path
from subprocess import check_output, CalledProcessError, STDOUT
import re, os

class Psi4Run:
    def __init__(self, inp, base='/dev/shm', verb=False):
        save_out = False #or True
        self.verb = verb #or True
        self.inp = inp
        #base = "/qscratch/%s" % os.environ['USER']
        if save_out:
            self.tmpdir = None
            self.dirname = Path( mkdtemp(dir=base) )
        else:
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
            if self.verb:
                print(self.out)
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
        state = 0 # read any optimized geometry found
        opt = []
        with open(self.dirname / "output.dat", encoding='utf-8') as f:
            for line in f:
                if state == 0 and 'Final optimized geometry' in line:
                    opt = []
                    state = 1
                elif state == 1 and line.strip()[:8] == 'Geometry':
                    state = 2 # then keyword
                elif state == 2 and line.strip() == '':
                    state = 3 # then blank line
                elif state == 3:
                    if line.strip() == '':
                        state = 0
                        continue
                    tok = line.split()
                    if len(tok) != 4:
                        print("Syntax error reading optimized geom!")
                        opt = []
                        state = 0
                        continue
                    opt.append([tok[0], [float(t) for t in tok[1:]]])

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
            if len(opt) > 0:
                attr['optimized geometry'] = opt
        return attr

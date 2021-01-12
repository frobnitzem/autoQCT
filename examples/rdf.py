# This example computes the order statistics
# of the distance to nearest water hydrogen atom.
# It assumes that every water is a separate molecule
# with atom ordering: ['O', 'H', 'H'].
#
import numpy as np

dr = 0.05
bins = int(8.0/dr)+1
max_N = 5

class Data:
    def __init__(self):
        self.frames = 0
        self.counts = np.zeros((max_N, bins), np.float)
        self.n = 0

    # Add an unsorted list of distances to the respective (sorted) histograms
    def add(self, r):
        self.frames += 1
        if len(r) > self.n:
            self.n = len(r)
        r.sort()
        for i,d in enumerate(r):
            j = int(d/dr)
            if j < bins:
                self.counts[i,j] += 1

    def cdf(self):
        N = self.frames + (self.frames == 0)
        g = self.counts[:self.n].cumsum(1) / N # cumulative dist'n
        return (np.arange(bins)+1)*dr, g

    def rdf(self):
        N = self.frames + (self.frames == 0)
        g = self.counts[:self.n] / N # radial dist'n
        V = (4*np.pi/3.0) * (np.arange(bins+1)*dr)**3
        g /= V[1:]-V[:-1]
        return (np.arange(bins)+0.5)*dr, g

x0 = Data()

dist = lambda x,y: np.sum((x-y)**2, -1)**0.5
def fn(qct, x, u):
    rad = []
    sys = qct.system
    k = sys[0].atoms()
    cm = np.sum(x[:k],0) / float(k)
    for m in sys[1:]:
        n = m.atoms()
        if n == 3 and ''.join(m.names) == 'OHH':
            rad.append(min( dist(cm, x[k+1])
                          , dist(cm, x[k+2])
                          ))
        k += n
    u.add(rad)

def output(name, u):
    print("%s has %d frames and %d max waters"%(name, u.frames, u.n))
    fname = "%s-rdf.xvg"%(name.rsplit('.',1)[0])
    r, M = u.rdf()
    with open(fname, "w") as f:
        for i in range(bins):
            f.write(" ".join( ["%g"%r[i]]
                            + ["%g"%y for y in M[:,i]] )
                    + "\n")

    fname = "%s-cdf.xvg"%(name.rsplit('.',1)[0])
    r, M = u.cdf()
    with open(fname, "w") as f:
        for i in range(bins):
            f.write(" ".join( ["%g"%r[i]]
                            + ["%g"%y for y in M[:,i]] )
                    + "\n")


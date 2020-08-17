from .psi4run import Psi4Run
from .molecule import *

# Input template
rrho_in = """#! rrho

molecule mol {{
{sys}
units angstrom
}}

set {{
  basis {basis}
}}

optimize('{theory}')
ks_e = frequencies('{theory}', dertype=1)
# TODO: check positivity of ks_e[6:]
"""

# Example output section:
"""
  Total E, Electronic energy at  298.15 [K]      -76.39670697 [Eh]
  Total H, Enthalpy at  298.15 [K]               -76.39576279 [Eh]
  Total G, Free enthalpy at  298.15 [K]          -76.41721220 [Eh]

~> gets formatted to {
     'H': -76.39670697,
     'E': -76.39576279,
     'G': -76.41721220
   }
"""
expr = r'\s*Total ([EHG]),.* 298.15 \[K\]\s*(.*) \[Eh\]'

def run_rrho(sys, crds, theory, basis):
    p = Psi4Run( rrho_in.format(sys=sys.fmt_psi4(crds, False),
                                theory=theory, basis=basis) )
    p.run()
    if p.err is None:
        print(p.scrape(expr))
    else:
        print("Error running psi4")

def test_run():
    sys = Sys([ Mol(['O', 'H', 'H'], 0, 1) ])
    crds = [[ 0.000, 0.000, 0.000],
            [ 0.803, 0.000, 0.596],
            [-0.803, 0.000, 0.596]]
    run_rrho(sys, crds, theory='wB97X-D', basis='aug-cc-pvdz')

#test_run()

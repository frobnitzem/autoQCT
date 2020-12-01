from .psi4run import Psi4Run
from .molecule import *

# EM Input template
em_in = """#! emin

molecule mol {{
{sys}
units angstrom
}}

set {{
  basis {basis}
}}

optimize('{theory}')
"""

# RRHO Input template
rrho_in = em_in + """
ks_e = frequencies('{theory}', dertype=1)
# TODO: check positivity of ks_e[6:]
"""

# Example emin output section:
"""
        Writing optimization data to binary file.
        Final energy is   -613.1645497661005
        Final (previous) structure:
        Cartesian Geometry (in Angstrom)
           CL    -0.7227130530  -1.0855534697  -0.5100229342
            O     0.1517684748   0.9073835467   1.9213631669
            H     0.7724497940   1.2272580628   1.2534531198
            H    -0.2985037558   0.2055096943   1.4156811211
            O     1.3087316090   1.2338330262  -0.8528720711
            H     0.7926017743   0.3925288301  -0.8076266338
            H     0.6304800221   1.8578656461  -1.1227950546
        Saving final (previous) structure.
"""

# Example rrho output section:
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

def run_rrho(sys, crds, theory, basis):
    expr = r'\s*Total ([EHG]),.* 298.15 \[K\]\s*(.*) \[Eh\]'
    p = Psi4Run( rrho_in.format(sys=sys.fmt_psi4(crds, False),
                                theory=theory, basis=basis) )
    p.run()

    if p.err is None:
        ret = p.scrape(expr, True)
    else:
        print("Error running psi4")
        ret = None
    return ret

def run_emin(sys, crds, theory, basis):
    expr = r'\s*Final (energy) is\s*(.*)'
    p = Psi4Run( em_in.format(sys=sys.fmt_psi4(crds, False),
                              theory=theory, basis=basis) )
    p.run()

    if p.err is None:
        ret = p.scrape(expr, True)
    else:
        print("Error running psi4")
        ret = None
    return ret['energy']

def test_run():
    sys = Sys([ Mol(['O', 'H', 'H'], 0, 1) ])
    crds = [[ 0.000, 0.000, 0.000],
            [ 0.803, 0.000, 0.596],
            [-0.803, 0.000, 0.596]]
    run_rrho(sys, crds, theory='wB97X-D', basis='aug-cc-pvdz')

#test_run()

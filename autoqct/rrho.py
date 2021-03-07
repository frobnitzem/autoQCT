from .psi4run import Psi4Run
from .molecule import *
import numpy as np

# EM Input template
em_in = """#! emin

molecule mol {{
{sys}
units angstrom
}}

set {{
  basis {basis}
  OPT_COORDINATES CARTESIAN
  geom_maxiter 100
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
    if len(crds) == 1:
        return run_rrho1(sys, crds, theory, basis)
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

def run_rrho1(sys, crds, theory, basis):
    E = run_emin(sys, crds, theory, basis, False)

    # STP
    T = 298.15
    P = 1.01325e5 # Pa (= 1 atm = default for psi4)
    # convert units!
    Na = 6.02214076e23 # pt / mol
    M = elem_dict[sys[0].names[0]][1]*1e-3/Na # kg
    kT = 8.3145*T / 6.02214076e23 # J
    h  = 6.62607015e-34 # J s
    # energy scale
    RT = 8.3145e-3*T/4.184 # kcal/mol

    # kg * J / J^2 s^2 = kg/s^2 / (kg m^2 / s^2) = 1/m^2
    # J/Pa = N*m / (N/m^2) = m^3
    Qt = (2*np.pi*M*kT/h**2)**1.5 * (kT/P)
    return {
             'H': E + 2.5*RT,
             'E': E,
             'G': E - RT*np.log(Qt),
             'optimized geometry': [
                [ sys[0].names[0], [ 0.0,0.0,0.0] ]
             ]
           }

# Also runs single-point energy calculations if optimize==False
def run_emin(sys, crds, theory, basis, optimize=True):
    expr = r'\s*Final (energy) is\s*(.*)'
    if optimize:
        tpl = em_in
    else:
        tpl = em_in.replace("optimize", "energy")
    p = Psi4Run( tpl.format(sys=sys.fmt_psi4(crds, False),
                            theory=theory, basis=basis) )
    p.run()

    if p.err is not None:
        print("Error running psi4")
        return None

    if not optimize:
        expr = '\s*Total (Energy) =\s*(.*)'
        ret = p.scrape(expr, True)
        return ret['Energy']

    ret = p.scrape(expr, True)
    return ret['energy']

def test_run(argv):
    sys = Sys([ Mol(['O', 'H', 'H'], 0, 1) ])
    crds = [[ 0.000, 0.000, 0.000],
            [ 0.803, 0.000, 0.596],
            [-0.803, 0.000, 0.596]]
    ans = run_rrho(sys, crds, theory=argv[1], basis=argv[2])
    print("Water:")
    print(ans)
    return

    for ion in ['F', 'Cl', 'Br', 'I']:
        sys = Sys([ Mol([ion], -1, 1) ])
        en = run_emin(sys, crds[:1], theory=argv[1], basis=argv[2],
                      optimize=False)
        print("%s-: %f"%(ion, en))

if __name__=="__main__":
    import sys
    test_run(sys.argv)

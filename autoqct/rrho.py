from psi4run import Psi4Run

# Input template
rrho_in = """#! rrho

molecule mol {{
{charge} {spin}
symmetry c1
{crds}
units angstrom
no_reorient
no_com
}}

set {{
  basis {basis}
}}

optimize('{theory}')
ks_e = frequencies('{theory}', dertype=1)
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

def run_rrho(crds, charge, spin, theory, basis):
    p = Psi4Run(rrho_in.format(crds=crds, charge=charge,
                               spin=spin, theory=theory,
                               basis=basis))
    p.run()
    if p.err is None:
        print(p.scrape(expr))
    else:
        print("Error running psi4")

def test_run():
    crds = "\n".join([
                "O  0.000 0.000 0.000",
                "H  0.803 0.000 0.596",
                "H -0.803 0.000 0.596"])
    run_rrho(crds, charge=0, spin=1, theory='b3lyp', basis="cc-pVDZ")

#test_run()

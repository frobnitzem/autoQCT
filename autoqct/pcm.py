from psi4run import Psi4Run

# Input template
pcm_in = """#! pcm

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
  scf_type pk
  pcm true
  pcm_scf_type total
}}

pcm = {{
    Units = Angstrom
    Medium {{
        SolverType = IEFPCM
        Solvent = Water
    }}

    Cavity {{
        RadiiSet = UFF
        Type = GePol
        Scaling = false
        Area = 0.3
        Mode = Implicit
    }}
}}

energy = energy('{theory}')
"""

# Example output section:
"""
  @RKS Final Energy:   -76.42663379788245

   => Energetics <=

    Nuclear Repulsion Energy =              8.7962296684390751
    One-Electron Energy =                -122.5818166091975741
    Two-Electron Energy =                  44.9184051361068342
    DFT Exchange-Correlation Energy =      -7.5501770580040963
    Empirical Dispersion Energy =           0.0000000000000000
    VV10 Nonlocal Energy =                  0.0000000000000000
    PCM Polarization Energy =              -0.0092749352266872
    Total Energy =                        -76.4266337978824453

~> gets formatted to {
     'Nuclear Repulsion': 8.7962296684390751,
     'One-Electron': -122.5818166091975741,
     ...
     'PCM Polarization': -0.0092749352266872,
     'Total Energy': -76.4266337978824453
   }
"""

def run_pcm(crds, charge, spin, theory, basis):
    p = Psi4Run(pcm_in.format(crds=crds, charge=charge,
                              spin=spin, theory=theory, basis=basis))
    p.run()
    if p.err is None:
        print(p.scrape(r'\s*(.*)\sEnergy\s*=\s*(.*)'))
    else:
        print("Error running psi4")

def test_run():
    crds = "\n".join([
                "O  0.000 0.000 0.000",
                "H  0.803 0.000 0.596",
                "H -0.803 0.000 0.596"])
    run_pcm(crds, charge=0, spin=1, theory='b3lyp', basis="cc-pVDZ")

#test_run()

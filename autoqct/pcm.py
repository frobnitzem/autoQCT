from .psi4run import Psi4Run
from .molecule import *

# Input template
pcm_in = """#! pcm

molecule mol {{
{sys}
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
        # how to manually set radii?
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

def run_pcm(sys, crds, theory, basis):
    p = Psi4Run(pcm_in.format(sys=sys.fmt_psi4(crds, False),
                              theory=theory, basis=basis))
    p.run()
    ret = None
    if p.err is None:
        ret = p.scrape(r'\s*(.*)\sEnergy\s*=\s*(.*)')
        #print(ret)
    else:
        print("Error running psi4")
    return ret

def test_run():
    sys = Sys([ Mol(['O', 'H', 'H'], 0, 1) ])
    crds = [[ 0.000, 0.000, 0.000],
            [ 0.803, 0.000, 0.596],
            [-0.803, 0.000, 0.596]]
    ret = run_pcm(sys, crds, theory='wB97X-D', basis='aug-cc-pvdz')
    print(ret)

#test_run()

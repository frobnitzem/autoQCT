from multiprocessing import Process, Pipe
from tempfile import TemporaryDirectory
import os
from .molecule import *

# FIXME: read the # of cores and memory from cmd-line args / config somewhere
# return counterpoise-corrected binding energy in kcal/mol
def calcE(conn, sys1, crd1, sys2, crd2, theory, basis):
    import psi4
    olddir = os.getcwd()
    tmpdir = TemporaryDirectory()
    os.chdir(tmpdir.name)
    psi4.core.set_output_file('output.dat', False)
    psi4.set_memory('100 GB')
    #psi4.set_num_threads(18)
    psi4.set_num_threads(1)

    psi4.set_options({'freeze_core': 'true'})
    cplx = psi4.geometry("{sys1}\n--\n{sys2}\nunits angstrom".format(
                  sys1=sys1.fmt_psi4(crd1), sys2=sys2.fmt_psi4(crd2))
           )
    de = psi4.energy('%s/%s'%(theory,basis), bsse_type='cp', molecule=cplx)
    conn.send( de * psi4.constants.hartree2kcalmol )
    os.chdir(olddir)
    conn.close()

def run_ebind(sys1, crd1, sys2, crd2, theory, basis):
    recv, send = Pipe(False)
    p = Process(target=calcE, args=(send,
                sys1, crd1, sys2, crd2, theory, basis))
    p.start()
    de = recv.recv()
    p.join() # this blocks until the process terminates
    print("Delta E = %g kcal/mol"%de)
    return de

def test_run():
    sys = Sys([ Mol(['O', 'H', 'H'], 0, 1) ])
    crd1 = [[ 0.000, 0.000, 0.000],
            [ 0.803, 0.000, 0.596],
            [-0.803, 0.000, 0.596]]
    crd2 = [[ 0.000, 3.000, 0.000],
            [ 0.803, 3.000, 0.596],
            [-0.803, 3.000, 0.596]]
    # or theory = ccsd(t), etc. or basis = aug-cc-pvdz
    run_ebind(sys, crd1, sys, crd2, theory='wB97X-D', basis='aug-cc-pvdz')

#test_run()

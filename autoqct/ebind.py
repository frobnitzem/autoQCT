from multiprocessing import Process, Pipe

def my_function(x):
    return x+100

# return counterpoise-corrected binding energy in kcal/mol
def calcE(conn, chg1, spin1, geom1, chg2, spin2, geom2, theory, basis):
    import psi4
    psi4.set_options({'freeze_core': 'true'})
    cplx = psi4.geometry("""{chg1} {spin1}
{geom1}
--
{chg2} {spin2}
{geom2}
units angstrom
""")
    de = psi4.energy('ccsd(t)/aug-cc-pvdz', bsse_type='cp', molecule=cplx)
    conn.send( de * psi4.constants.hartree2kcalmol )
    conn.close()

def run_ebind(charge1, spin1, crds1, charge2, spin2, crds2, theory, basis):
    recv, send = Pipe(False)
    p = Process(target=calcE, args=(send,
		chg1, spin1, crds1, chg2, spin2, crds2, theory, basis))
    p.start()
    ret = recv.recv()
    ret = p.join() # this blocks until the process terminates
    return ret

def test_run():
    crds1 = "\n".join([
                "O  0.000 0.000 0.000",
                "H  0.803 0.000 0.596",
                "H -0.803 0.000 0.596"])
    crds2 = "\n".join([
                "O  0.000 3.000  0.000",
                "H  0.803 3.000 -0.596",
                "H -0.803 3.000 -0.596"])
    # or theory = ccsd(t), etc.
    run_ebind(0, 1, crds1, 0, 1, crds2, theory='b3lyp', basis='aug-cc-pvdz')

#test_run()

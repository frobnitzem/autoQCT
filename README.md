# autoQCT <!--![Build Status][build] ![Documentation Status][docs]-->

This project aims to automate most steps in the calculation
of solvation free energies using quasi-chemical theory.

# Quasi-Chemical Theory

QCT treats each ion-plus-n-water cluster as a molecular species of the system under analysis, then provides a concise format,

![equation](https://latex.codecogs.com/gif.latex?\mu_{X^-}^{(ex)}&space;=&space;-&space;RT&space;\ln&space;K_n^{(0)}&space;&plus;&space;\rho_{\mathrm{H}_2\mathrm{O}}^{~~~~n}&space;&plus;&space;RT&space;\ln&space;p_{X^-}(n)&space;&plus;&space;\left\(\mu_{(\mathrm{H}_2\mathrm{O})_n}^{(ex)}&space;-&space;n&space;\mu_{\mathrm{H}_2\mathrm{O}}^{(ex)}&space;\right&space;\))

for the solvation free energy of the ion, X-, based on a Born-Haber cycle for solvating its n-water hydrate.
This requires estimating the solvation free energy of the hydrate, as well
as the probability that the probability that X- is surrounded by n waters
in solution.  The upshot is that that the chemical contribution of
forming the cluster itself can be computed using high-level quantum
mechanical calculations in gas phase.

## Prerequisites

This package deals specifically with the parallel
calculation of cluster energetics needed to compute the
![equation](https://latex.codecogs.com/gif.latex?K_n^{(0)}) term.

Before starting, you must have a CP2K trajectory
of an inner-shell solute/solvent cluster.
You are free to define what exactly is meant
(spatially) by a cluster
geometry that possesses the "inner-shell" property.

We define "an inner-shell cluster
as one where all ligands of X are inner-shell", and, further,
"a water is inner-shell if it has an H atom within a distance
of λ from any atom in X."
This is a natural answer for anions, where the variety
of H-bond donation structures that can form with nearby water
would have complicated using hydrogen-bond metrics.

We suggest determining the "cut-off" distance, λ,
for interacting H atoms through radial distributions
provided from *ab initio* molecular dynamics in bulk.
This package does not do that step at present.

## Isolated Cluster Energetics

To evaluate cluster energetics, and to address the discrepancy in harmonic approximation and free energy experiments, 
we select configurations containing n water molecules within λ of the ion of interest.

These configurations go through a molecule deletion scheme to evaluate the free energy,

![equation](https://latex.codecogs.com/gif.latex?\Delta&space;U&space;=&space;E(\gamma_n\sigma)&space;-&space;E(\gamma_{n-1}\sigma)&space;-&space;E(\gamma&space;\sigma)&space;&plus;&space;E(\sigma))

so that the formation energy of the n-water hydrate can be determined
by sequential water additions,

![equation](https://latex.codecogs.com/gif.latex?n&space;K_n^{(0)}&space;=&space;\frac{K_1^{(0)}K_{n-1}^{(0)}}{\langle&space;e^{\beta\Delta&space;U}\rangle}_n)

Where γ denote water ligands and σ is the ion of interest.
This naturally requires an estimate of K(n-1).  That quantity
can be calculated by using the rigid-rotor harmonic oscillator (RRHO)
approximation or from water addition to K(n-2).
Plotting these two estimates as a function of n shows when the
RRHO approximation has converged.

## Operation

This set of scripts is run from your own python script that
imports the autoqct package.  You describe your
calculation's input parameters (like your soute/solvent's
definition of inner-shell-ness) using a class
that extends 'QCT'.  Running autoQCT on that class
will work like `make` to create a dictory structure
that manages the calculation's state.

Like make, the autoQCT function can be run with a target
so that it will stop at a particular stage of the process.
It can also be set to verbose or quiet operation mode.

Internally, the state proceeds through three steps:

  1. Validation of the input trajectory (valid target)

  2. Job-script creation (script target)

    This will create (frames)*(n+1) sub-directories: one for each
    trajectory frame - with all n waters and all-but-one water.
    It will also create an example launch script to
    submit all jobs *en bloc*.

    It will also create four sub-directories to launch
    energy minimization jobs for water, ion, n-water ion,
    and n-1-water ion complexes.

  3. Output validation (output target)

    This step will validate all outputs from the job.
    On detecting errors, it will create a partial launch script
    to re-submit all incomplete jobs.

  4. Output collection and summary (Kn target)

    This step will compute the necessary averages to report a final estimate of
    ![equation](https://latex.codecogs.com/gif.latex?n&space;\ln&space;K_n^{(0)}).

# How do I use it?

First, verify that you have the following dependencies installed
(in order of increasing difficulty):

* Linux-like environment with bash
* python 3 with numpy
* [CP2K](https://www.cp2k.org)
* [Psi4](http://www.psicode.org/)
  * A method that works well is to use the cmake build process in
    a python virtualenv.

Second, download and install the project:

```bash
git clone https://github.com/frobnitzem/autoQCT.git
cd autoQCT
python3 setup.py install
```

Third, obtain a CP2K trajectory (`.xyz`) of the
solute/solvent cluster you are interested in. 
Put it in an empty directory and call it something like `cp2k_cluster.xyz`.

Fourth, create a `qct.py` input file following the examples in
the [examples directory][2].
Test out your cluster-ness filter by running
```
./qct.py test_filter cp2k_cluster.xyz 0
```
This will run your `clustered` function on frame 0 and show
the coordinates and the filter result.  This way you can check
whether your filter is working as expected.

When you are satisfied with the filter, it's time
to create the dataset.  Do this by executing:
```
./qct.py init cp2k_cluster.xyz
```

If all goes well, you should now see lots of `frame_nnn.npy` files containing
coordinates of your cluster.  You should also have a `targets.yaml`
file listing out QM calculations that need to be run.

Submit those to the batch cluster using as many processors as you
have available.  To do this, copy both `scripts/batch.sh` and
and `scripts/targets.yaml` to the current directory.  Adapt
`scripts/batch.sh` to your use case, and run it.  You should be able
to check job progress from its log file.  Re-run this until
it says all jobs are completed.

At this point, you should have all the following output files:

* `dg_nm1_n.en`
* `dg_solv.txt`
* `nm1.yaml`


Finally, if this work has been useful in creating a published work,
please acknowledge it by citing [Quasi-chemical theory for anion hydration and specific ion effects: Cl-(aq) vs. F-(aq), Chem. Phys. Lett. X. 4:100037, 2019.][1].

# Contributing

We welcome contributions to this project via pull-requests or
as [issues/bug reports](https://github.com/frobnitzem/autoQCT/issues/new).

[1]: https://doi.org/10.1016/j.cpletx.2019.100037
[2]: https://github.com/frobnitzem/autoQCT/blob/master/examples
[3]: http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library

<!--
[build]: https://travis-ci.org/frobnitzem/autoqct.svg?branch=develop
[docs]: https://readthedocs.org/projects/autoqct/badge/?version=latest
-->


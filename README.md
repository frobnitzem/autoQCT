# autoQCT <!--![Build Status][build] ![Documentation Status][docs]-->

This project aims to automate most steps in the calculation
of solvation free energies using quasi-chemical theory.

# How do I use it?

First, verify that you have the following dependencies installed
(in order of increasing diffulty):

* Linux-like environment with bash
* python 3 with numpy
* [Gromacs XTC library][3]
  - `git clone git://git.gromacs.org/libxdrfile.git`
  - `cd libxdrfile && ./configure && make install`
  - optionally use `--prefix` option to configure to customize install location
* [CP2K](https://www.cp2k.org)
* [Gaussian 09 or 16](https://gaussian.com) with ADMP dynamics

Second, download and install the project:

```bash
git clone https://github.com/frobnitzem/autoQCT.git
cd autoQCT
python3 setup.py install
```

Third, obtain a gromacs trajectory (`.trr` or `.xtc`) of the
solute/solvent system you are interested in. 

Fourth, create an input file following the examples in
the [examples directory][2], and execute it with `python3 my-file.py`.
You should see status output as the calculation progresses.

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


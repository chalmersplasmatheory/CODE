[![Documentation Status](https://readthedocs.org/projects/code/badge/?version=latest)](https://code.readthedocs.io/en/latest/?badge=latest)

# CODE -- Collisional Distributions of Electrons v2
This repository contains the new, object-oriented implementation of the 0D-2P
kinetic solver CODE (Collision Distributions of Electrons). The original solver
CODE has been developed gradually over several years in the plasma theory group
at Chalmers University of Technology, Sweden, in order to simulate the momentum
space dynamics of runaway electrons. This version of the code is a refactoring
of the original code, utilising the object orientation features of Matlab.

This version of the code was written by Albert Johansson.

## References
If publishing results computed with CODE, please cite the following papers:

[1] Landreman M., Stahl A. and Fülöp T., [*Comp. Phys. Comm.* **185** (3) 847-855](https://doi.org/10.1016/j.cpc.2013.12.004) (2013).

[2] Stahl A., Embréus O., Papp G., Landreman M. and Fülöp T., [*Nucl. Fusion* **56** (11) 112009](https://doi.org/10.1088/0029-5515/56/11/112009) (2016).

BibTeX:
```
@article {Landreman2013,
    author = {M. Landreman and A. Stahl and T. Fülöp},
    title = {Numerical calculation of the runaway electron distribution function and associated synchrotron emission},
    journal = {Computer Physics Communications},
    volume = {185},
    number = {3},
    pages = {847 - 855},
    year = {2014},
    issn = {0010-4655},
    doi = {10.1016/j.cpc.2013.12.004},
    url = {https://doi.org/10.1016/j.cpc.2013.12.004}
}
@article {Stahl2016,
	title = {Kinetic modelling of runaway electrons in dynamic scenarios},
	author = {A. Stahl and O. Embr{\'{e}}us and G. Papp and M. Landreman and T. Fülöp},
	year = 2016,
	month = {jul},
	publisher = {{IOP} Publishing},
	volume = {56},
	number = {11},
	pages = {112009},
	journal = {Nuclear Fusion},
	doi = {10.1088/0029-5515/56/11/112009},
	url = {https://doi.org/10.1088/0029-5515/56/11/112009}
}
```

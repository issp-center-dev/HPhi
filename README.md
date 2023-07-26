HΦ
====

A numerical solver package for a wide range of quantum lattice models including Hubbard-type itinerant electron hamiltonians, quantum spin models, and Kondo-type hamiltonians for itinerant electrons coupled with quantum spins. The Lanczos algorithm for finding ground states and newly developed Lanczos-based algorithm for finite-temperature properties of these models are implemented for parallel computing (hybrid parallelization with OpenMP and MPI). A broad spectrum of users including experimental researchers is cordially welcome.

### Methods
Lanczos algorithm, thermal pure quantum state, full diagonalization  

### Target models
Hubbard model, Heisenberg model, Kondo lattice model, Kitaev model, Kitaev-Heisenberg model, multi-orbital Hubbard model

### Available physical quantities
specific heat, susceptibility, ground state energy, free energy, structure factors


## Requirement
C compiler (intel, Fujitsu, GNU, etc. )  
LAPACK library (intel MKL, Fujitsu, ATLAS, etc.)  
MPI library (If you do not use MPI, it is not necessary.)

## Install

You can install HΦ and also get a manual for HΦ from a [release note](https://github.com/issp-center-dev/HPhi/releases).

## Pre-Installed system
- System B of The Institute Solid State Physics jpint-use supercomputer (cleck [here](http://www.issp.u-tokyo.ac.jp/supercom/visitor/x92nxz/hphi) for details).

## Licence

The distribution of the program package and the source codes for HPhi follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).We hope that you cite the following references when you publish the results using HΦ (hphi):

[“Quantum lattice model solver HΦ”, M. Kawamura, K. Yoshimi, T. Misawa, Y. Yamaji, S. Todo, and N. Kawashima, Computer Physics Communications 217, 180 (2017).](https://doi.org/10.1016/j.cpc.2017.04.006)

[“Update of HΦ: Newly added functions and methods in versions 2 and 3”, K. Ido, M. Kawamura, Y. Motoyama, K. Yoshimi, Y. Yamaji, S. Todo, N. Kawashima, and T. Misawa, arXiv:2307.13222.](https://arxiv.org/abs/2307.13222)

Bibtex:

@article{KAWAMURA2017180,
title = {Quantum lattice model solver HΦ},
journal = {Computer Physics Communications},
volume = {217},
pages = {180-192},
year = {2017},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2017.04.006},
url = {https://www.sciencedirect.com/science/article/pii/S0010465517301200},
author = {Mitsuaki Kawamura and Kazuyoshi Yoshimi and Takahiro Misawa and Youhei Yamaji and Synge Todo and Naoki Kawashima}
}

@misc{ido2023update,
      title={Update of $\mathcal{H}\Phi$: Newly added functions and methods in versions 2 and 3},
      author={Kota Ido and Mitsuaki Kawamura and Yuichi Motoyama and Kazuyoshi Yoshimi and Youhei Yamaji and Synge Todo and Naoki Kawashima and Takahiro Misawa},
      year={2023},
      eprint={2307.13222},
      archivePrefix={arXiv},
      primaryClass={cond-mat.str-el}
}


## Official page
- [HPhi portal site](https://www.pasums.issp.u-tokyo.ac.jp/hphi/en)
- [User-Manual](https://www.pasums.issp.u-tokyo.ac.jp/hphi/en/doc/manual)
- [HPhi tutorial](https://issp-center-dev.github.io/HPhi/manual/develop/tutorial/en/html/index.html)
- [HPhi Gallery](https://isspns-gitlab.issp.u-tokyo.ac.jp/hphi-dev/hphi-gallery)

## Author
Youhei Yamaji, Takahiro Misawa, Synge Todo, Kota Ido, Yuichi Motoyama, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Naoki Kawashima.

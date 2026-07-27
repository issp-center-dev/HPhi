<p align="center">
  <img src="doc/figs/logo_HPhi.png" width="300">
</p>

<h1 align="center">HΦ</h1>

<p align="center">
  <a href="https://github.com/issp-center-dev/HPhi/releases"><img src="https://img.shields.io/github/v/release/issp-center-dev/HPhi" alt="GitHub release"></a>
  <a href="https://github.com/issp-center-dev/HPhi/actions/workflows/main.yml"><img src="https://github.com/issp-center-dev/HPhi/actions/workflows/main.yml/badge.svg" alt="Build Status"></a>
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="License: GPL v3"></a>
  <a href="https://doi.org/10.1016/j.cpc.2017.04.006"><img src="https://img.shields.io/badge/DOI-10.1016%2Fj.cpc.2017.04.006-blue" alt="DOI"></a>
  <a href="https://www.pasums.issp.u-tokyo.ac.jp/hphi/en/doc/manual"><img src="https://img.shields.io/badge/docs-manual-brightgreen" alt="Manual"></a>
  <a href="https://issp-center-dev.github.io/HPhi/manual/develop/api/html/index.html"><img src="https://img.shields.io/badge/docs-API-orange" alt="API Reference"></a>
</p>

A numerical solver package for a wide range of quantum lattice models including Hubbard-type itinerant electron hamiltonians, quantum spin models, and Kondo-type hamiltonians for itinerant electrons coupled with quantum spins.

## Features

- **High-performance parallel computing**: Hybrid parallelization with OpenMP and MPI
- **Multiple solvers**: Lanczos algorithm, LOBPCG, thermal pure quantum (TPQ) state, full diagonalization
- **Wide range of models**: Hubbard, Heisenberg, Kondo, Kitaev, multi-orbital systems
- **Rich output**: Ground state energy, free energy, specific heat, susceptibility, structure factors

## Quick Start

### Requirements

- C compiler (Intel, GCC, Fujitsu, etc.)
- LAPACK library (Intel MKL, OpenBLAS, etc.)
- MPI library (optional, for parallel computation)
- Optional high-performance backends for full diagonalization (ScaLAPACK, MAGMA, ELPA): see the installation guide in the user manual

### Installation

```bash
# Clone the repository
git clone https://github.com/issp-center-dev/HPhi.git
cd HPhi

# Build with CMake
mkdir build && cd build
cmake ..
make -j4

# Run tests
ctest
```

For detailed installation instructions, see the [release notes](https://github.com/issp-center-dev/HPhi/releases).

### Basic Usage

```bash
# Generate input files using StdFace
HPhi -s input.in

# Or run directly with expert mode
HPhi -e namelist.def
```

## Supported Models

| Model Type | Examples |
|------------|----------|
| Hubbard | Single-band, multi-orbital Hubbard models |
| Spin | Heisenberg, XY, Ising, Kitaev, Kitaev-Heisenberg |
| Kondo | Kondo lattice model |
| t-J | t-J model |
| Spinless Fermion | Spinless fermion models |

## Methods

| Method | Description |
|--------|-------------|
| **Lanczos** | Ground state and low-lying excited states |
| **LOBPCG** | Block Lanczos for multiple eigenvalues |
| **TPQ** | Finite-temperature properties via thermal pure quantum states |
| **Full Diagonalization** | All eigenvalues for small systems |
| **Real-time Evolution** | Time-dependent simulations |

## Documentation

- [User Manual (English)](https://www.pasums.issp.u-tokyo.ac.jp/hphi/en/doc/manual)
- [User Manual (Japanese)](https://www.pasums.issp.u-tokyo.ac.jp/hphi/ja/doc/manual)
- [Tutorial](https://issp-center-dev.github.io/HPhi/manual/develop/tutorial/en/html/index.html)
- [API Reference (Doxygen)](https://issp-center-dev.github.io/HPhi/manual/develop/api/html/index.html)
- [HPhi Gallery](https://isspns-gitlab.issp.u-tokyo.ac.jp/hphi-dev/hphi-gallery)

## Pre-installed Systems

- System B at The Institute for Solid State Physics (ISSP), The University of Tokyo ([details](http://www.issp.u-tokyo.ac.jp/supercom/visitor/x92nxz/hphi))

## Citation

If you use HΦ in your research, please cite:

**Version 1:**
> M. Kawamura, K. Yoshimi, T. Misawa, Y. Yamaji, S. Todo, and N. Kawashima,
> "Quantum lattice model solver HΦ",
> [Computer Physics Communications **217**, 180 (2017)](https://doi.org/10.1016/j.cpc.2017.04.006)

**Versions 2 & 3:**
> K. Ido, M. Kawamura, Y. Motoyama, K. Yoshimi, Y. Yamaji, S. Todo, N. Kawashima, and T. Misawa,
> "Update of HΦ: Newly added functions and methods in versions 2 and 3",
> [Computer Physics Communications **298**, 109093 (2024)](https://doi.org/10.1016/j.cpc.2024.109093)

<details>
<summary>BibTeX</summary>

```bibtex
@article{KAWAMURA2017180,
  title   = {Quantum lattice model solver {H$\Phi$}},
  author  = {Mitsuaki Kawamura and Kazuyoshi Yoshimi and Takahiro Misawa and
             Youhei Yamaji and Synge Todo and Naoki Kawashima},
  journal = {Computer Physics Communications},
  volume  = {217},
  pages   = {180--192},
  year    = {2017},
  doi     = {10.1016/j.cpc.2017.04.006}
}

@article{IDO2024109093,
  title   = {Update of {$\mathcal{H}\Phi$}: Newly added functions and methods in versions 2 and 3},
  author  = {Kota Ido and Mitsuaki Kawamura and Yuichi Motoyama and
             Kazuyoshi Yoshimi and Youhei Yamaji and Synge Todo and
             Naoki Kawashima and Takahiro Misawa},
  journal = {Computer Physics Communications},
  volume  = {298},
  pages   = {109093},
  year    = {2024},
  doi     = {10.1016/j.cpc.2024.109093}
}
```

</details>

## License

HΦ is distributed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

## Authors

Youhei Yamaji, Takahiro Misawa, Synge Todo, Kota Ido, Yuichi Motoyama, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Naoki Kawashima

## Links

- [HPhi Portal Site](https://www.pasums.issp.u-tokyo.ac.jp/hphi/en)
- [GitHub Repository](https://github.com/issp-center-dev/HPhi)
- [Issue Tracker](https://github.com/issp-center-dev/HPhi/issues)

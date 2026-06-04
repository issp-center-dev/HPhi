# Off-diagonal dynamical Green function: reference-value validation

These small exact-diagonalization (ED) scripts independently reproduce the
reference values hard-coded in the off-diagonal dynamical Green-function
regression tests, so the references can be re-derived and audited rather than
trusted blindly.

The tests compute `G_BA(z) = <gs| B^dagger (z-H)^-1 A |gs>` (ket excitation `A`,
bra excitation `B`) via HPhi's shifted BiCG; these scripts compute the same
quantity by dense diagonalization.

## Scripts

| script | model | validates tests |
|---|---|---|
| `exact_diag_hubbard_offdiag.py` | 4-site Hubbard ring (t=1, U=4, half filling) | `spectrum_hubbard_offdiag_single`, `spectrum_hubbard_offdiag_pair_density`, `spectrum_hubbard_offdiag_single_mpi` |
| `exact_diag_spin_offdiag.py`    | 4-site Heisenberg ring (J=1), spin S=1/2 and S=1 | `spectrum_spin_offdiag_szsz`, `spectrum_spingc_offdiag_szsz` |

## Requirements / how to run

Only `numpy` is required. Activate the project Python environment first:

```sh
pyenv activate venv312
python test/tools/exact_diag_hubbard_offdiag.py
python test/tools/exact_diag_spin_offdiag.py
```

Each script prints the ground-state energy `E0` (which equals the `OmegaOrg`
used in the test inputs) followed by `omega  Re G  Im G` for `omega = 0.0 .. 0.4`.
Compare those columns to the `reference.dat` blocks embedded in the matching
`test/*.sh`.

## Agreement and tolerances

- HPhi uses shifted BiCG with `LanczosEps = 14`, so its output matches the exact
  ED here to about `1e-8` (single / Sz-Sz) and `1e-6` (pair-density). The serial
  regression tests embed the HPhi output itself, so their internal `L2 diff`
  is `0.000000`; the ED scripts are the *independent* correctness check.
- The MPI test (`spectrum_hubbard_offdiag_single_mpi`) compares an `np=4` run to
  the serial/exact reference with a `< 1e-5` tolerance, because MPI changes the
  BiCG convergence path (observed `L2 diff ~ 3e-7`).

## Conventions (must match HPhi)

- `z = omega + OmegaOrg + i*OmegaIm`, with `OmegaOrg = E0`.
- `numpy.vdot(B, x)` conjugates its first argument, matching HPhi's
  `VecProdMPI(<B phi|, .)` projection, i.e. the `<B phi|` direction.
- The bra operators are written in the `*Bra` files as `B` (not `B^dagger`); the
  computed quantity is `<phi| B^dagger (z-H)^-1 A |phi>`.
- `Sz_i` is the standard pair-operator representation
  (S=1/2: `-0.5 n_{i,0} + 0.5 n_{i,1}`; S=1: `-n_{i,0} + n_{i,2}`). The overall
  sign convention is irrelevant for the `Sz`-`Sz` correlation.

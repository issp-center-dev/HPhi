"""Side-effect-free canonical sector basis shared by sector test suites.

HPhi's phase convention normalizes sum_g conj(chi(g)) T_g |r> for the
smallest raw state r of each orbit, with a positive real coefficient at r,
and orders the columns by r.
"""
import numpy as np


def canonical_basis(projector):
    """Canonical U from an exact sector projector, one column per orbit."""
    covered = np.zeros(len(projector), dtype=bool)
    columns = []
    for r in range(len(projector)):
        if covered[r] or projector[r, r].real < 1e-12:
            continue
        v = projector[:, r].copy()
        v /= np.linalg.norm(v)
        # Identity for an exact Hermitian projector; it fixes the convention
        # against rounding in the diagonal element.
        v *= np.conj(v[r]) / abs(v[r])
        covered |= np.abs(v) > 1e-12
        columns.append(v)
    if not columns:
        return np.empty((len(projector), 0), dtype=complex)
    return np.column_stack(columns)


def check_basis(projector, basis):
    """Assert orthonormality, completeness, phase and column order of U."""
    np.testing.assert_allclose(basis.conj().T @ basis,
                               np.eye(basis.shape[1]), rtol=0, atol=1e-12)
    if projector is not None:
        np.testing.assert_allclose(basis @ basis.conj().T, projector,
                                   rtol=0, atol=1e-12)
    representatives = []
    for column in basis.T:
        r = np.flatnonzero(abs(column) > 1e-12)[0]
        assert column[r].real > 0 and abs(column[r].imag) < 1e-12, column[r]
        representatives.append(r)
    assert all(a < b for a, b in zip(representatives, representatives[1:])), \
        representatives

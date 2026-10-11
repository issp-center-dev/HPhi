"""Requested sector FullDiag correlations for every eigenstate, all formats."""
from pathlib import Path
import shutil
import sys

import numpy as np

import symmetry_fulldiag_common as common

ATOL = 1e-8


def expected_values(case, vectors):
    return {family: np.array([[np.vdot(v, oq @ v) for _, oq in requests]
                              for v in vectors.T])
            for family, requests in case.operators.items()}


def compare(actual, expected, label):
    np.testing.assert_allclose(actual.real, expected.real, rtol=0, atol=ATOL, err_msg=label)
    np.testing.assert_allclose(actual.imag, expected.imag, rtol=0, atol=ATOL, err_msg=label)


def read_all(case, aggregate):
    output = case.path/'output'
    return {family: common.read_green(output, family, len(case.hq), aggregate,
                                      [row for row, _ in requests])
            for family, requests in case.operators.items()}


def cluster_means(case, values):
    return [values[start:end].mean(axis=0) for start, end in common.cluster_edges(case.hq)]


def reference_cluster_means(case):
    """Basis-independent cluster average: trace(P_c O)/dim_c."""
    _, reference = np.linalg.eigh(case.hq)
    expected = expected_values(case, reference)
    return {family: cluster_means(case, values) for family, values in expected.items()}


def run(case, hphi, solver, launcher, label, vectors, aggregate):
    if (case.path/'output').exists():
        shutil.rmtree(case.path/'output')
    flags = {'OutputEigenVec': int(vectors), 'OutputGreenFormat': int(aggregate)}
    text = common.run_fd(case, hphi, solver, flags, launcher, label)
    return text, case.path/'output'


def green_files(output):
    return [file for pattern in ('*cisajs*', '*ThreeBody*', '*FourBody*', '*SixBody*', '*NBodyG*')
            for file in output.glob(pattern)]


def check_case(case, hphi, solver, label, launcher, ranks, with_correlations, plain_sym):
    dim = len(case.hq)
    tag = f'{case.name} {label} np={ranks}'
    # Correlations and states together: each row against the saved state.
    _, output = run(case, hphi, solver, launcher, f'both_{label}', True, False)
    energies = np.loadtxt(output/'zvo_energy_sector.dat', ndmin=2)[:, 1]
    vectors, _ = common.read_states(output, dim, ranks)
    common.check_eigensystem(case.hq, energies, vectors)
    expected = expected_values(case, vectors)
    actual = read_all(case, False)
    for family in case.operators:
        compare(actual[family], expected[family], f'{tag} {family} normal')
    references = reference_cluster_means(case)
    for family, values in actual.items():
        for mean, reference in zip(cluster_means(case, values), references[family]):
            compare(mean, reference, f'{tag} {family} cluster mean')
    manifest = common.read_manifest(output)
    assert manifest['completed_outputs'] == 'eigenvalues,correlations,eigenvectors', manifest
    assert manifest['fulldiag_status'] == 'complete'
    # Aggregate rows of the same run kind against that run's own states.
    _, output = run(case, hphi, solver, launcher, f'aggregate_{label}', True, True)
    vectors, _ = common.read_states(output, dim, ranks)
    expected = expected_values(case, vectors)
    aggregate = read_all(case, True)
    for family in case.operators:
        compare(aggregate[family], expected[family], f'{tag} {family} aggregate')
        assert not list(output.glob(f'zvo_{common.GREEN_STEMS[family]}_eigen0.dat'))
    # Correlations only: degenerate states may differ between runs.
    _, output = run(case, hphi, solver, launcher, f'correlations_{label}', False, False)
    assert not list(output.glob('*eigenvec*'))
    only = read_all(case, False)
    for family, values in only.items():
        for mean, reference in zip(cluster_means(case, values), references[family]):
            compare(mean, reference, f'{tag} {family} correlations-only cluster mean')
    manifest = common.read_manifest(output)
    assert manifest['completed_outputs'] == 'eigenvalues,correlations'
    # States only and neither: no Green output at all.
    (case.path/'sym.def').write_text(plain_sym)
    _, output = run(case, hphi, solver, launcher, f'vectors_{label}', True, False)
    assert not green_files(output) and len(list(output.glob('*eigenvec*'))) == dim*ranks
    _, output = run(case, hphi, solver, launcher, f'none_{label}', False, False)
    assert not green_files(output) and not list(output.glob('*eigenvec*'))
    assert common.read_manifest(output)['output_scope'] == 'eigenvalues'
    (case.path/'sym.def').write_text(with_correlations)
    clusters = len(common.cluster_edges(case.hq))
    print(f'{tag}: {dim} states x {sum(len(r) for r in case.operators.values())} rows, '
          f'{clusters} clusters, normal/aggregate/only/none PASS', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
    assert mode == 'all', mode
    launcher, ranks = common.launcher_from_env()
    selection = common.selected_cases()
    root = common.suite_root(f'correlations-{selection}')
    backends = common.solver_backends(launcher, ranks)
    cases = common.cases(root, selection)
    if selection == 'all':
        cases.append(common.zero_hamiltonian_case(root))
    exercised = {}
    for case in cases:
        plain_sym = (case.path/'sym.def').read_text()
        if not common.attach_correlations(case, False):
            print(f'{case.name}: no independent operator oracle; correlations skipped', flush=True)
            continue
        with_correlations = (case.path/'sym.def').read_text()
        for solver, label, backend_launcher in backends:
            backend_ranks = common.spingc.mpi_size(backend_launcher)
            if solver == 3 and len(case.hq) < common.grid_max(backend_ranks):
                continue  # the grid rejection is covered by the states suite
            check_case(case, hphi, solver, label, backend_launcher, backend_ranks,
                       with_correlations, plain_sym)
            exercised.setdefault(label, []).append((case.name, len(case.hq), backend_ranks))
    models = {case.model for case in cases}
    for label, runs in exercised.items():
        print(f'{label} exercised: ' + ', '.join(f'{name}(D={dim},np={n})' for name, dim, n in runs),
              flush=True)
    print(f'sector FullDiag correlations ({selection}): {len(models)} models PASS', flush=True)


if __name__ == '__main__':
    main()

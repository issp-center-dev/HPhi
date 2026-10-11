"""Sector FullDiag I/O failures must stop every rank without a hang."""
import os
from pathlib import Path
import shutil
import sys

import symmetry_fulldiag_common as common

TIMEOUT = 30


def blocked_run(case, hphi, solver, launcher, label, blocked, message, flags):
    """Pre-create an output name as a directory so that only its opener fails."""
    output = case.path/'output'
    if output.exists():
        shutil.rmtree(output)
    output.mkdir()
    (output/blocked).mkdir()
    text = common.run_fd(case, hphi, solver, flags, launcher, label, message, timeout=TIMEOUT)
    manifest = (output/'zvo_symmetry_sector.dat').read_text()
    assert 'fulldiag_status=complete' not in manifest, manifest
    return text


def green_mode(hphi, launcher, ranks):
    root = common.suite_root('failures-green')
    cases = {case.name: case for case in common.cases(root, 'all')}
    # Hubbard requests every Green family, so later TwoBodyG files exist.
    case = cases['Hubbard_L4_k0']
    assert common.attach_correlations(case, False)
    families = set(case.operators)
    assert {'one', 'two', 'three', 'four', 'six', 'nbody'} <= families, families
    scenarios = [
        ('onebody', 'zvo_cisajs_eigen0.dat', 0),
        ('twobody_later_file', 'zvo_ThreeBody_eigen0.dat', 0),
        ('twobody_last_file', 'zvo_SixBody_eigen0.dat', 0),
        ('nbody', 'zvo_NBodyG_eigen0.dat', 0),
        ('second_state', 'zvo_cisajs_eigen1.dat', 0),
        ('aggregate', 'zvo_cisajs_eigen.dat', 1),
        ('aggregate_twobody', 'zvo_cisajscktalt_eigen.dat', 1),
    ]
    backends = common.solver_backends(launcher, ranks)
    if ranks > 1 and len(backends) == 1:
        print('green failures not exercised: a multi-rank FullDiag needs ScaLAPACK or ELPA, '
              'which this build lacks', flush=True)
        common.require_backends({})
    for solver, label, backend_launcher in backends:
        if ranks > 1 and solver == 0:
            continue  # LAPACK is serial; the collective cases need several ranks
        for name, blocked, aggregate in scenarios:
            blocked_run(case, hphi, solver, backend_launcher, f'{name}_{label}', blocked,
                        'sector correlation output open failed',
                        {'OutputEigenVec': 0, 'OutputGreenFormat': aggregate})
            print(f'{label} np={common.spingc.mpi_size(backend_launcher)} {name}: '
                  f'{blocked} rejected on every rank', flush=True)


# PR-A stages: (stage, rank, backends, diagnostic, flags). Allocation and
# checkpoint stages fail on rank 1, root-only output stages on rank 0.
IO_STAGES = (
    ('vector_alloc', 1, ('ScaLAPACK', 'ELPA'),
     'sector FullDiag state buffer allocation failed', {'OutputEigenVec': 1}),
    ('matrix_alloc', 1, ('ELPA',),
     'sector Hamiltonian column buffer allocation failed', {'OutputEigenVec': 1}),
    ('checkpoint_write', 1, ('ScaLAPACK', 'ELPA'),
     'symmetry checkpoint write failed', {'OutputEigenVec': 1}),
    ('checkpoint_close', 1, ('ScaLAPACK', 'ELPA'),
     'symmetry checkpoint write failed', {'OutputEigenVec': 1}),
    ('checkpoint_rename', 1, ('ScaLAPACK', 'ELPA'),
     'symmetry checkpoint publish failed', {'OutputEigenVec': 1}),
    ('green_close', 0, ('ScaLAPACK', 'ELPA'),
     'sector correlation output close failed', {'OutputEigenVec': 0}),
)


def io_mode(hphi, launcher, ranks):
    """Injected rank-local I/O and allocation failures end every rank.

    matrix_alloc is the ELPA column-panel buffer; it needs an ELPA build."""
    if ranks < 2:
        print('io failures not exercised: injection on rank 1 needs an MPI launcher with '
              'at least two ranks', flush=True)
        common.require_backends({})  # a job that requires backends must not skip
        return
    backends = {label: solver for solver, label, _ in common.solver_backends(launcher, ranks)}
    if set(backends) == {'LAPACK'}:
        print('io failures not exercised: a multi-rank FullDiag needs ScaLAPACK or ELPA, '
              'which this build lacks', flush=True)
        common.require_backends({})
        return
    probe = common.probe_executable(hphi)
    root = common.suite_root('failures-io')
    cases = {case.name: case for case in common.cases(root, 'all')}
    case = cases['Hubbard_L4_k0']
    assert common.attach_correlations(case, False)
    exercised = []
    for stage, rank, labels, message, flags in IO_STAGES:
        usable = [label for label in labels if label in backends]
        if not usable:
            print(f'{stage} not exercised: needs {"/".join(labels)} in this build', flush=True)
            continue
        label = usable[-1]
        env = {'HPHI_TEST_FD_FAIL': stage, 'HPHI_TEST_FD_FAIL_RANK': str(rank)}
        output = case.path/'output'
        if output.exists():
            shutil.rmtree(output)
        common.run_fd(case, probe, backends[label], {**flags, 'OutputGreenFormat': 0},
                      launcher, f'{stage}_{label}', message, timeout=TIMEOUT, env_extra=env)
        manifest = (output/'zvo_symmetry_sector.dat').read_text()
        assert 'fulldiag_status=complete' not in manifest, (stage, manifest)
        if stage.startswith('checkpoint'):
            # The failed state is never published, and no temporary remains.
            assert not (output/f'zvo_eigenvec_0_rank_{rank}.dat').exists(), stage
            assert not list(output.glob('*.part')), stage
        exercised.append((stage, label))
        print(f'{stage} on rank {rank} with {label} np={ranks}: every rank stopped', flush=True)
        # The same environment does not reach production HPhi.
        if output.exists():
            shutil.rmtree(output)
        common.run_fd(case, hphi, backends[label], {**flags, 'OutputGreenFormat': 0},
                      launcher, f'{stage}_{label}_production', timeout=TIMEOUT, env_extra=env)
        assert common.read_manifest(output)['fulldiag_status'] == 'complete'
    assert exercised, 'no I/O failure stage was exercised'
    if os.environ.get('HPHI_HAS_ELPA') == '1':
        assert ('matrix_alloc', 'ELPA') in exercised


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'green'
    launcher, ranks = common.launcher_from_env()
    if mode == 'green':
        green_mode(hphi, launcher, ranks)
    elif mode == 'io':
        io_mode(hphi, launcher, ranks)
    else:
        raise ValueError(mode)
    print(f'sector FullDiag failures ({mode}) PASS', flush=True)


if __name__ == '__main__':
    main()

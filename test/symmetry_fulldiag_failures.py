"""Sector FullDiag I/O failures must stop every rank without a hang."""
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
    assert {'one', 'two', 'three', 'nbody'} <= families, families
    scenarios = [
        ('onebody', 'zvo_cisajs_eigen0.dat', 0),
        ('twobody_later_file', 'zvo_ThreeBody_eigen0.dat', 0),
        ('nbody', 'zvo_NBodyG_eigen0.dat', 0),
        ('second_state', 'zvo_cisajs_eigen1.dat', 0),
        ('aggregate', 'zvo_cisajs_eigen.dat', 1),
        ('aggregate_twobody', 'zvo_cisajscktalt_eigen.dat', 1),
    ]
    backends = common.solver_backends(launcher, ranks)
    for solver, label, backend_launcher in backends:
        if ranks > 1 and solver == 0:
            continue  # LAPACK is serial; the collective cases need several ranks
        for name, blocked, aggregate in scenarios:
            blocked_run(case, hphi, solver, backend_launcher, f'{name}_{label}', blocked,
                        'sector correlation output open failed',
                        {'OutputEigenVec': 0, 'OutputGreenFormat': aggregate})
            print(f'{label} np={common.spingc.mpi_size(backend_launcher)} {name}: '
                  f'{blocked} rejected on every rank', flush=True)


def main():
    hphi = Path(sys.argv[1]).resolve()
    mode = sys.argv[2] if len(sys.argv) > 2 else 'green'
    launcher, ranks = common.launcher_from_env()
    if mode == 'green':
        green_mode(hphi, launcher, ranks)
    else:
        raise ValueError(mode)
    print(f'sector FullDiag failures ({mode}) PASS', flush=True)


if __name__ == '__main__':
    main()

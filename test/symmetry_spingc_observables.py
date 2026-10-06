"""SpinGC symmetry-sector magnetization moments and iterative solvers."""
import os
from pathlib import Path
import re
import shlex
import sys
import tempfile

import numpy as np
import symmetry_spingc_common as c


def translation(nsite, momentum):
    permutations = [[(site+shift) % nsite for site in range(nsite)]
                    for shift in range(nsite)]
    characters = np.exp(-2j*np.pi*momentum*np.arange(nsite)/nsite)
    return permutations, characters


def fixture(label):
    if label == 'A':
        nsite = 8
        permutations, characters = translation(nsite, 0)
        families = {'Trans': [(i, a, i, b, .5, 0)
                              for i in range(nsite)
                              for a, b in [(1, 0), (0, 1)]]}
        hamiltonian = -sum(c.spin_operators(nsite)['x'])
    else:
        nsite = 6 if label == 'D' else 8
        momentum = 1 if label in ('B1', 'D') else 0
        permutations, characters = translation(nsite, momentum)
        families = c.mixed_families(nsite)
        hamiltonian = c.mixed_hamiltonian(nsite)
    basis, representatives, _ = c.sector_basis(nsite, permutations, characters)
    sector_hamiltonian = basis.conj().T @ hamiltonian @ basis
    return (nsite, permutations, characters, families, basis,
            representatives, sector_hamiltonian)


def rank_slice(dimension, ranks, rank):
    count, remainder = divmod(dimension, ranks)
    size = count + (rank < remainder)
    offset = count*rank + min(rank, remainder)
    return offset, size


def write_probe_inputs(path, vector, ranks):
    for rank in range(ranks):
        offset, size = rank_slice(len(vector), ranks, rank)
        text = ''.join('{:.17g} {:.17g}\n'.format(value.real, value.imag)
                       for value in vector[offset:offset+size])
        (path/'probe-input.rank{}.dat'.format(rank)).write_text(text)


def expected_moments(vector, representatives, nsite):
    magnetization = np.array([bin(int(rep)).count('1')-nsite/2
                              for rep in representatives])
    return np.array([np.vdot(vector, magnetization*vector).real,
                     np.vdot(vector, magnetization*magnetization*vector).real])


def parse_moments(text):
    sz = re.search(r'^SpinGCProbe Sz ([^\n]+)$', text, re.M)
    sz2 = re.search(r'^SpinGCProbe Sz2 ([^\n]+)$', text, re.M)
    assert sz and sz2, text
    result = np.array([float(sz.group(1)), float(sz2.group(1))])
    assert np.isfinite(result).all()
    return result


def probe_moments(path, probe, launcher, layout, vector, representatives,
                  nsite, expected_error=None, extra_env=None):
    ranks = c.mpi_size(launcher)
    write_probe_inputs(path, vector, ranks)
    env = dict(os.environ, HPHI_TEST_SPINGC_ACTION='moments',
               HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    if extra_env:
        env.update(extra_env)
    text = c.run_case(path, probe, 'moments', launcher, env, expected_error)
    if expected_error is None:
        np.testing.assert_allclose(parse_moments(text),
                                   expected_moments(vector, representatives, nsite),
                                   atol=1e-8, rtol=0)
    return text


def parse_energy(path):
    text = (path/'output/zvo_energy.dat').read_text()
    match = re.search(r'^\s*Energy\s+([^\s]+)', text, re.M)
    assert match, text
    return float(match.group(1)), text


def check_eigenvector(vector, energy, hamiltonian, representatives, nsite):
    np.testing.assert_allclose(np.vdot(vector, vector), 1, atol=1e-9, rtol=0)
    residual = np.linalg.norm(hamiltonian@vector-energy*vector)
    scale = max(1, np.linalg.norm(hamiltonian, 2)*np.linalg.norm(vector))
    assert residual/scale <= 1e-8, (residual, scale)
    return expected_moments(vector, representatives, nsite)


def checkpoint_vector(path, dimension):
    files = list((path/'output').glob('zvo_eigenvec_0_rank_*.dat'))
    assert files
    return c.join_rank_vectors(files, dimension)


def final_lanczos_vector(path, dimension, representatives, ranks):
    vector = np.empty(dimension, complex)
    covered = np.zeros(dimension, bool)
    files = list(path.glob('spingc-lanczos.rank*.dat'))
    assert len(files) == ranks
    for file in files:
        for line in file.read_text().splitlines():
            index, representative, real, imag = line.split()
            index = int(index)-1
            assert 0 <= index < dimension and not covered[index]
            assert int(representative) == int(representatives[index])
            vector[index] = complex(float(real), float(imag))
            covered[index] = True
    assert covered.all() and np.isfinite(vector).all()
    return vector


def arbitrary_vector(dimension, seed):
    index = np.arange(1, dimension+1, dtype=float)
    vector = np.cos(index*(.17+seed*.01)) + 1j*np.sin(index*(.29+seed*.02))
    return vector/np.linalg.norm(vector)


def exercise_fixture(root, label, hphi, probe, launcher, layout):
    (nsite, permutations, characters, families, basis, representatives,
     hamiltonian) = fixture(label)
    dimension = len(representatives)
    ranks = c.mpi_size(launcher)
    if label == 'D':
        assert dimension == 9
        assert sum(rank_slice(dimension, ranks, rank)[1] == 0
                   for rank in range(ranks)) == max(0, ranks-dimension)

    moments_path = root/'{}_{}_moments'.format(label, layout)
    c.write_case(moments_path, nsite, permutations, characters, families, 3, {})
    vector = arbitrary_vector(dimension, len(label))
    probe_moments(moments_path, probe, launcher, layout, vector,
                  representatives, nsite)
    if label == 'A':
        raw_ground = np.full(1 << nsite, 1/np.sqrt(1 << nsite), complex)
        analytic = basis.conj().T@raw_ground
        analytic /= np.linalg.norm(analytic)
        probe_moments(moments_path, probe, launcher, layout, analytic,
                      representatives, nsite)
        np.testing.assert_allclose(expected_moments(analytic, representatives, nsite),
                                   [0, 2], atol=1e-12, rtol=0)

    cg_path = root/'{}_{}_cg'.format(label, layout)
    c.write_case(cg_path, nsite, permutations, characters, families, 3,
                 {'CalcMod': {'OutputEigenVec': 1},
                  'ModPara': {'LanczosEps': 18}})
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    c.run_case(cg_path, hphi, 'cg', launcher, env)
    energy, energy_text = parse_energy(cg_path)
    vector = checkpoint_vector(cg_path, dimension)
    moments = check_eigenvector(vector, energy, hamiltonian,
                                representatives, nsite)
    output_sz = re.search(r'^\s*Sz\s+([^\s]+)', energy_text, re.M)
    assert output_sz
    np.testing.assert_allclose(float(output_sz.group(1)), moments[0], atol=1e-8, rtol=0)
    probe_moments(cg_path, probe, launcher, layout, vector,
                  representatives, nsite)
    if label == 'A':
        np.testing.assert_allclose([energy, moments[0], moments[1]],
                                   [-4, 0, 2], atol=3e-8, rtol=0)

    lanczos_path = root/'{}_{}_lanczos'.format(label, layout)
    c.write_case(lanczos_path, nsite, permutations, characters, families, 0,
                 {'CalcMod': {'CalcEigenVec': 1},
                  'ModPara': {'LanczosEps': 18}})
    lanczos_env = dict(env, HPHI_TEST_SPINGC_ACTION='lanczos')
    c.run_case(lanczos_path, probe, 'lanczos', launcher, lanczos_env)
    energy, energy_text = parse_energy(lanczos_path)
    vector = final_lanczos_vector(lanczos_path, dimension,
                                  representatives, ranks)
    moments = check_eigenvector(vector, energy, hamiltonian,
                                representatives, nsite)
    output_sz = re.search(r'^Sz\s+([^\s]+)', energy_text, re.M)
    assert output_sz
    np.testing.assert_allclose(float(output_sz.group(1)), moments[0], atol=1e-8, rtol=0)
    assert not list((lanczos_path/'output').glob('*eigenvec*'))
    print('{} {} dim={} CG/Lanczos/moments passed'.format(
        label, layout, dimension), flush=True)


def failures(root, probe, launcher, layout):
    nsite, permutations, characters, families, _, representatives, _ = fixture('D')
    dimension = len(representatives)
    vector = arbitrary_vector(dimension, 11)

    def case(name):
        path = root/'failure_{}_{}'.format(name, layout)
        c.write_case(path, nsite, permutations, characters, families, 3, {})
        write_probe_inputs(path, vector, c.mpi_size(launcher))
        return path

    for name, mutate in [
            ('short', lambda path: (path/'probe-input.rank0.dat').write_text('')),
            ('long', lambda path: (path/'probe-input.rank0.dat').write_text(
                (path/'probe-input.rank0.dat').read_text()+'0 0\n')),
            ('missing', lambda path: (path/'probe-input.rank0.dat').unlink()),
            ('nan', lambda path: (path/'probe-input.rank0.dat').write_text(
                'nan 0\n'+(path/'probe-input.rank0.dat').read_text().split('\n', 1)[1]))]:
        path = case(name)
        mutate(path)
        c.run_case(path, probe, name, launcher,
                   dict(os.environ, HPHI_TEST_SPINGC_ACTION='moments',
                        HPHI_SYMMETRY_BASIS_LAYOUT=layout),
                   'SpinGC moments probe failed')
    path = case('storage')
    victim = c.mpi_size(launcher)-1
    c.run_case(path, probe, 'storage', launcher,
               dict(os.environ, HPHI_TEST_SPINGC_ACTION='moments',
                    HPHI_SYMMETRY_BASIS_LAYOUT=layout,
                    HPHI_TEST_SPINGC_INVALID_STORAGE_RANK=str(victim)),
               'SpinGC moments probe failed')
    if c.mpi_size(launcher) == 1 and layout == 'replicated':
        path = root/'lanczos_output_gate'
        c.write_case(path, nsite, permutations, characters, families, 0,
                     {'CalcMod': {'CalcEigenVec': 1, 'OutputEigenVec': 1}})
        c.run_case(path, probe, 'output_gate', [], dict(os.environ),
                   'does not support OutputEigenVec with Lanczos')


def raw_regression(root, hphi):
    nsite, permutations, characters, families, _, _, _ = fixture('A')
    path = root/'raw_spingc'
    c.write_case(path, nsite, permutations, characters, families, 3,
                 {'ModPara': {'LanczosEps': 18}})
    lines = [line for line in (path/'sym.def').read_text().splitlines()
             if not line.startswith('TransSym ')]
    (path/'sym.def').write_text('\n'.join(lines)+'\n')
    text = c.run_case(path, hphi, 'raw', [], dict(os.environ))
    energy, values = parse_energy(path)
    sz = float(re.search(r'^\s*Sz\s+([^\s]+)', values, re.M).group(1))
    np.testing.assert_allclose([energy, sz], [-4, 0], atol=3e-8, rtol=0)
    assert 'TransSym' not in text


def main():
    hphi, probe = [Path(arg).resolve() for arg in sys.argv[1:]]
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = c.mpi_size(launcher)
    assert ranks in (1, 4, 16), 'observables suite requires explicit np1, np4, or np16'
    root = Path(tempfile.mkdtemp(prefix='symmetry_spingc_observables_', dir='.'))
    print('artifacts: {}'.format(root.resolve()), flush=True)
    for layout in ('replicated', 'distributed'):
        for label in ('A', 'B0', 'B1', 'D'):
            exercise_fixture(root, label, hphi, probe, launcher, layout)
        failures(root, probe, launcher, layout)
    if ranks == 1:
        raw_regression(root, hphi)


if __name__ == '__main__':
    main()

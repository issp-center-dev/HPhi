"""Independent tensor references and portable I/O for SpinGC sector tests."""
import os
from pathlib import Path
import struct
import subprocess

import numpy as np

E10 = np.array([[0, 0], [1, 0]], dtype=complex)
E01 = E10.conj().T
Sx, Sy = (E10 + E01)/2, (E10 - E01)/(2j)
Sz = np.diag([-.5, .5])


def tensor_at(nsite, site, matrix):
    result = np.ones((1, 1), complex)
    for pos in reversed(range(nsite)):
        result = np.kron(result, matrix if pos == site else np.eye(2))
    return result


def spin_operators(nsite):
    return {name: [tensor_at(nsite, i, matrix) for i in range(nsite)]
            for name, matrix in dict(plus=E10, minus=E01, x=Sx, y=Sy, z=Sz).items()}


def mixed_hamiltonian(nsite, chirality=.07):
    s = spin_operators(nsite)
    p, m, x, y, z = (s[key] for key in ('plus', 'minus', 'x', 'y', 'z'))
    h = np.zeros((2**nsite, 2**nsite), complex)
    for i in range(nsite):
        j = (i+1) % nsite
        h += (-.83*x[i] - .19*y[i] - .11*z[i] + .37*z[i]@z[j]
              + .23*(p[i]@m[j] + m[i]@p[j])
              + .17*(p[i]@p[j] + m[i]@m[j])
              + 1j*chirality*(p[i]@m[j] - m[i]@p[j])
              + .13*(x[i]@z[j] + z[i]@x[j]))
    return h


def sector_basis(nsite, permutations, characters):
    dim = 1 << nsite
    transforms = [[sum(((r >> i) & 1) << perm[i] for i in range(nsite))
                   for r in range(dim)] for perm in permutations]
    projector = np.zeros((dim, dim), complex)
    for moved, chi in zip(transforms, characters):
        projector[moved, np.arange(dim)] += np.conj(chi)/len(transforms)
    reps, columns = [], []
    for r in range(dim):
        if min(moved[r] for moved in transforms) != r:
            continue
        v = projector[:, r].copy()
        norm = np.linalg.norm(v)
        if norm < 1e-12:
            continue
        v /= norm
        v *= np.conj(v[r])/abs(v[r])
        reps.append(r)
        columns.append(v)
    basis = np.column_stack(columns) if columns else np.empty((dim, 0), complex)
    return basis, np.array(reps, dtype=np.uint64), projector


def interall_pair(indices, value):
    a, s, b, t, c, u, d, v = indices
    z = complex(value)
    return [tuple(indices) + (z.real, z.imag),
            (d, v, c, u, b, t, a, s, z.real, -z.imag)]


def mixed_interall_bond(i, j, chirality=.07):
    rows = interall_pair((i, 1, i, 0, j, 0, j, 1), 1j*chirality)
    for sigma, sign in ((1, 1), (0, -1)):
        rows += interall_pair((i, 1, i, 0, j, sigma, j, sigma), .0325*sign)
        rows += interall_pair((i, sigma, i, sigma, j, 1, j, 0), .0325*sign)
    assert len(rows) == 10
    return rows


def mixed_families(nsite, chirality=.07, all_interall=False):
    field = -.83*Sx-.19*Sy-.11*Sz
    families = {'Trans': [(i, a, i, b, -field[a, b].real, -field[a, b].imag)
                          for i in range(nsite) for a in range(2) for b in range(2)]}
    rows = []
    for i in range(nsite):
        j = (i+1) % nsite
        rows += mixed_interall_bond(i, j, chirality)
        if all_interall:
            for a in range(2):
                for b in range(2):
                    rows.append((i, a, i, a, j, b, j, b,
                                 .37*Sz[a, a]*Sz[b, b], 0))
            rows += interall_pair((i, 1, i, 0, j, 0, j, 1), .23)
            rows += interall_pair((i, 1, i, 0, j, 1, j, 0), .17)
    if not all_interall:
        for key, value in [('Ising', .37), ('Exchange', .23), ('PairLift', .17)]:
            families[key] = [(i, (i+1) % nsite, value) for i in range(nsite)]
    else:
        # Represent the field as E_i^{ab} (E_j^{00}+E_j^{11}).
        families.pop('Trans')
        for i in range(nsite):
            j = (i+1) % nsite
            for sigma in range(2):
                for a in range(2):
                    rows.append((i, a, i, a, j, sigma, j, sigma, field[a, a].real, 0))
                rows += interall_pair((i, 1, i, 0, j, sigma, j, sigma), field[1, 0])
    assert len(rows) == nsite*(26 if all_interall else 10)
    families['InterAll'] = rows
    return families


def validate_interall(rows):
    totals = {}
    for row in rows:
        key = tuple(row[:8])
        totals[key] = totals.get(key, 0) + complex(*row[8:])
    for key, value in totals.items():
        if key[1] == key[3] and key[5] == key[7]:
            assert abs(value.imag) < 1e-14
            continue
        a, s, b, t, c, u, d, v = key
        assert abs(totals.get((d, v, c, u, b, t, a, s), 0)-value.conjugate()) < 1e-14


def families_matrix(nsite, families):
    """Tensor interpretation of expert rows, without HPhi raw apply helpers."""
    units = {(i, a, b): tensor_at(nsite, i, np.eye(2)[:, a:a+1] @ np.eye(2)[b:b+1])
             for i in range(nsite) for a in range(2) for b in range(2)}
    s = spin_operators(nsite)
    h = np.zeros((1 << nsite, 1 << nsite), complex)
    for name, rows in families.items():
        if name == 'InterAll':
            validate_interall(rows)
        for row in rows:
            if name == 'Trans':
                i, a, j, b, re, im = row
                assert i == j
                h -= complex(re, im)*units[i, a, b]
            elif name == 'InterAll':
                i, a, ii, b, j, c, jj, d, re, im = row
                assert i == ii and j == jj
                h += complex(re, im)*units[i, a, b]@units[j, c, d]
            else:
                i, j, value = row
                if name == 'Ising':
                    # Reader: Hund=-J/2, CoulombInter=-J/4. Spin diagonal
                    # Hund contributes -Hund to equal spins, Coulomb is constant.
                    expanded = value/2*(units[i, 0, 0]@units[j, 0, 0]
                                        + units[i, 1, 1]@units[j, 1, 1])-value/4*np.eye(1 << nsite)
                    np.testing.assert_allclose(expanded, value*s['z'][i]@s['z'][j], atol=1e-14)
                    h += expanded
                elif name == 'Exchange':
                    h += value*(s['plus'][i]@s['minus'][j]+s['minus'][i]@s['plus'][j])
                elif name == 'PairLift':
                    h += value*(s['plus'][i]@s['plus'][j]+s['minus'][i]@s['minus'][j])
                else:
                    raise AssertionError(name)
    return h


def definition(path, name, rows, count=None, key='NItems'):
    body = '\n'.join(' '.join(map(str, row)) for row in rows)
    (path/name).write_text('====\n{} {}\n====\n====\n====\n{}\n'.format(
        key, len(rows) if count is None else count, body))


def write_case(path, nsite, permutations, characters, families, method, options):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    calc = dict(CalcType=method, CalcModel=4, OutputMode=0, OutputDataHead=1)
    calc.update(options.get('CalcMod', {}))
    (path/'calc.def').write_text(''.join('{} {}\n'.format(k, v) for k, v in calc.items()))
    mod = dict(Nsite=nsite, Lanczos_max=400, initial_iv=-1, exct=1,
               LanczosEps=12, LanczosTarget=1, LargeValue=100, PreCG=0)
    mod.update(options.get('ModPara', {}))
    (path/'mod.def').write_text('====\nModel_Parameters 0\n====\n====\n====\n'
        'CDataFileHead zvo\nCParaFileHead zqp\n====\n'
        + ''.join('{} {}\n'.format(k, v) for k, v in mod.items()))
    definition(path, 'loc.def', [(i, 1) for i in range(nsite)], nsite)
    rows = [(g, complex(chi).real, complex(chi).imag) for g, chi in enumerate(characters)]
    rows += [(g, i, target, 1) for g, perm in enumerate(permutations) for i, target in enumerate(perm)]
    definition(path, 'group.def', rows, len(permutations), 'NQPTrans')
    names = ['CalcMod calc.def', 'ModPara mod.def', 'LocSpin loc.def', 'TransSym group.def']
    for keyword, values in families.items():
        if keyword == 'InterAll':
            validate_interall(values)
        definition(path, keyword+'.def', values)
        names.append('{} {}.def'.format(keyword, keyword))
    (path/'sym.def').write_text('\n'.join(names)+'\n')


def mpi_size(launcher):
    if not launcher:
        return 1
    for flag in ('-np', '-n', '--np', '--n'):
        if flag in launcher:
            return int(launcher[launcher.index(flag)+1])
    return None


def run_case(path, executable, label, launcher, env, expected_error=None):
    launcher = list(launcher)
    if expected_error is not None and mpi_size(launcher) == 1:
        launcher = []
    result = subprocess.run(launcher + [str(Path(executable).resolve()), '-e', 'sym.def'],
                            cwd=path, env=env, capture_output=True, text=True, timeout=120)
    text = result.stdout + result.stderr
    (Path(path)/(label+'.log')).write_text(text)
    if expected_error is None:
        assert result.returncode == 0, text
    else:
        assert result.returncode != 0 and expected_error in text, text
    return text


def read_checkpoint(path):
    data = Path(path).read_bytes()
    assert len(data) >= 28*8
    header = struct.unpack('<28Q', data[:28*8])
    assert header[:4] == (0x0a31565349485048, 1, 1, 128)
    assert len(data) == 28*8+16*header[15]
    vector = np.frombuffer(data[28*8:], dtype='<c16').copy()
    assert np.isfinite(vector).all()
    return header, vector


def join_rank_vectors(paths, dimension):
    parts = [read_checkpoint(p) for p in paths]
    parts.sort(key=lambda part: (part[0][14], part[0][12]))
    result = np.empty(dimension, complex)
    offset = 0
    ranks = set()
    for header, vector in parts:
        assert header[10] == dimension and header[11] == len(parts)
        assert header[12] not in ranks and header[12] < len(parts)
        ranks.add(header[12])
        assert header[14] == offset and offset+len(vector) <= dimension
        result[offset:offset+len(vector)] = vector
        offset += len(vector)
    assert offset == dimension and np.isfinite(result).all()
    return result


def checkpoint_digest(parts):
    """Independent version-1 payload digest, including empty rank ownership.

    parts contains (header, vector) pairs; hash binary64 bits, not rounded text.
    This deliberately does not validate norm/finite values, so rejection tests
    can isolate those guards with otherwise consistent checksums.
    """
    mask = (1 << 64)-1
    xor_hash, sum_hash = 0, 0
    for header, vector in parts:
        data = struct.pack('<3Q', header[12], header[14], header[15])
        data += np.asarray(vector, dtype='<c16').tobytes()
        value = 14695981039346656037
        for byte in data:
            value = ((value ^ byte)*1099511628211) & mask
        xor_hash ^= value
        sum_hash = (sum_hash+value) & mask
    return xor_hash, sum_hash


if __name__ == '__main__':
    raise SystemExit('Import this module from a SpinGC test suite.')

"""Portable expert input and probe I/O; the physics oracle lives in reference.py."""
import os
from pathlib import Path
import shlex
import numpy as np
from symmetry_spingc_common import definition, run_case, mpi_size
from symmetry_kondo_reference import Case


def merged_rows(rows: list[tuple[tuple[int, ...], complex]]) -> list[tuple[tuple[int, ...], complex]]:
    totals = {}
    for key, value in rows:
        totals[key] = totals.get(key, 0j) + value
    if not all(np.isfinite(z) for z in totals.values()):
        raise ValueError('nonfinite merged coefficient')
    return [(key, totals[key]) for key in sorted(totals) if totals[key] != 0]


def expert_case(path: Path, case: Case, *, method: int, options: dict,
                empty: bool = False) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    p = case.cells
    local = list(range(p)) if case.layout == 'block' else list(range(0, 2*p, 2))
    cond = list(range(p, 2*p)) if case.layout == 'block' else list(range(1, 2*p, 2))
    calc = dict(CalcType=method, CalcModel=5 if case.model == 'KondoGC' else 2,
                OutputMode=0, OutputDataHead=1)
    calc.update(options.get('CalcMod', {}))
    (path/'calc.def').write_text(''.join(f'{k} {v}\n' for k, v in calc.items()))
    mod = dict(Nsite=2*p, Lanczos_max=400, initial_iv=-1, exct=1,
               LanczosEps=12, LanczosTarget=1, LargeValue=100, PreCG=0)
    if case.ncond is not None:
        mod['Ncond'] = case.ncond
    if case.sz2 is not None:
        mod['2Sz'] = case.sz2
    mod.update(options.get('ModPara', {}))
    mod = {k: v for k, v in mod.items() if v is not None}
    (path/'mod.def').write_text('====\nModel_Parameters 0\n====\n====\n====\n'
        'CDataFileHead zvo\nCParaFileHead zqp\n====\n'
        + ''.join(f'{k} {v}\n' for k, v in mod.items()))
    definition(path, 'loc.def', [(i, int(i in local)) for i in range(2*p)], p)
    rows = [(g, np.cos(-2*np.pi*case.momentum*g/p), np.sin(-2*np.pi*case.momentum*g/p)) for g in range(p)]
    rows += [(g, sites[j], sites[(j+g) % p], 1) for g in range(p)
             for sites in (local, cond) for j in range(p)]
    definition(path, 'group.def', rows, p, 'NQPTrans')
    with (path/'group.def').open('a') as f:
        f.write(f'# MomentumIndex {case.momentum}\n')
    families = {}
    if not empty:
        transfer, inter = [], []
        def adjoint(key):
            pairs = list(zip(key[::2], key[1::2]))
            return tuple(x for pair in reversed(pairs) for x in pair)
        def pair(rows, key, value):
            rows.extend([(key, complex(value)), (adjoint(key), complex(value).conjugate())])
        for j, (l, c) in enumerate(zip(local, cond)):
            ll, cc = local[(j+1) % p], cond[(j+1) % p]
            for spin, sz in ((0, .5), (1, -.5)):
                # HPhi Transfer carries an overall minus sign.
                pair(transfer, (cc, spin, c, spin), .73*np.exp(.17j))
                transfer.extend([((l, spin, l, spin), .11*sz),
                                 ((c, spin, c, spin), .11*.7*sz+.23)])
                for other, sz_other in ((0, .5), (1, -.5)):
                    inter.extend([((c, spin, c, spin, cc, other, cc, other), .13),
                                  ((l, spin, l, spin, c, other, c, other), 1.10*sz*sz_other),
                                  ((l, spin, l, spin, ll, other, ll, other), .19*sz*sz_other)])
            inter.append(((c, 0, c, 0, c, 1, c, 1), .41))
            pair(inter, (l, 0, l, 1, c, 1, c, 0), .83/2)
            if case.model != 'Kondo':
                pair(transfer, (l, 0, l, 1), .27/2)
                pair(transfer, (c, 0, c, 1), -.08j)
                pair(inter, (l, 0, l, 1, ll, 0, ll, 1), .09)
        for key, values in [('Trans', transfer), ('InterAll', inter)]:
            families[key] = [(*indices, value.real, value.imag) for indices, value in merged_rows(values)]
    families.update(options.get('families', {}))
    names = ['CalcMod calc.def', 'ModPara mod.def', 'LocSpin loc.def', 'TransSym group.def']
    for key, rows in families.items():
        definition(path, key+'.def', rows)
        names.append(f'{key} {key}.def')
    (path/'sym.def').write_text('\n'.join(names)+'\n')
    return path


def read_manifest(path: Path) -> dict[str, str]:
    path = Path(path)
    if path.is_dir():
        path = path/'output/zvo_symmetry_sector.dat'
    return dict(line.split('=', 1) for line in path.read_text().splitlines())


def read_vector_parts(path: Path, stem: str, dimension: int) -> np.ndarray:
    vector = np.empty(dimension, complex)
    seen = np.zeros(dimension, bool)
    files = sorted(Path(path).glob(stem+'*.dat'))
    assert files, (path, stem)
    for file in files:
        for line in file.read_text().splitlines():
            index, rep, real, imag = line.split()
            index = int(index)-1
            assert 0 <= index < dimension and not seen[index]
            seen[index] = True
            vector[index] = complex(float(real), float(imag))
    assert seen.all() and np.isfinite(vector).all()
    return vector


def run_probe(path: Path, case: Case, executable: Path, *, action: str,
              layout: str, vector: np.ndarray | None = None) -> dict:
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks is not None
    if vector is not None:
        # Contiguous block ownership gives the first remainder ranks one extra row.
        for rank in range(ranks):
            base, remainder = divmod(len(vector), ranks)
            start = base*rank+min(rank, remainder)
            stop = start+base+(rank < remainder)
            np.savetxt(path/f'probe-input.rank{rank}.dat',
                       np.column_stack((vector[start:stop].real, vector[start:stop].imag)), fmt='%.17g')
    env = dict(os.environ, HPHI_TEST_SYMMETRY_ACTION=action, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    # A previous invocation must not supply missing output files.
    for stale in path.glob('sector_probe_rank_*'):
        stale.unlink()
    text = run_case(path, executable, action, launcher, env)
    data = read_manifest(path)
    dim = int(data['sector_dim'])
    files = sorted(path.glob('sector_probe_rank_*.dat'))
    assert len(files) == ranks
    infos = [read_manifest(file.with_suffix('.info')) for file in files]
    assert all(int(info['dim']) == dim for info in infos)
    assert sum(int(info['local_dim']) for info in infos) == dim
    for info in infos:
        assert int(info['raw_basis_elements']) == int(info['raw_diagonal_elements']) == 0
    result = dict(manifest=data, info=infos, text=text)
    if action == 'matvec':
        matrix, seen = np.empty((dim, dim), complex), np.zeros((dim, dim), bool)
        representatives = np.zeros(dim, dtype=np.uint64)
        for file in files:
            info = read_manifest(file.with_suffix('.info'))
            offset, size = int(info['offset']), int(info['local_dim'])
            lines = file.read_text().splitlines()
            assert len(lines) == dim*size
            for line in lines:
                row, col, rep, real, imag = line.split()
                row, col = int(row)-1, int(col)-1
                assert offset <= row < offset+size and 0 <= col < dim and not seen[row, col]
                seen[row, col] = True
                representatives[row] = int(rep)
                matrix[row, col] = complex(float(real), float(imag))
        assert seen.all() and np.isfinite(matrix).all()
        result.update(matrix=matrix, representatives=representatives)
    else:
        result['vector'] = read_vector_parts(path, 'sector_probe_rank_', dim)
        representatives = np.zeros(dim, dtype=np.uint64)
        for file in files:
            for line in file.read_text().splitlines():
                index, rep, *_ = line.split()
                representatives[int(index)-1] = int(rep)
        result['representatives'] = representatives
    return result


def read_kondo_checkpoint(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the Kondo-only v2 wire format without changing the v1 reader."""
    data = Path(path).read_bytes()
    assert len(data) >= 30*8, 'truncated Kondo checkpoint header'
    header = np.frombuffer(data[:30*8], dtype='<u8').copy()
    assert tuple(header[:4]) == (0x0a31565349485048, 2, 2, 128), tuple(header[:4])
    assert len(data) == 30*8+int(header[15])*16
    payload = np.frombuffer(data[30*8:], dtype='<c16').copy()
    assert np.isfinite(payload).all()
    return header, payload

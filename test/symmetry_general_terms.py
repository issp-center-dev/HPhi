"""Complete spectra against tensor-product operators, raw HPhi and signed projectors."""
import cmath
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import numpy as np

HPHI = str(Path(sys.argv[1]).resolve())
ROOT = Path("symmetry_general_terms")
ROOT.mkdir(exist_ok=True)
MPI = shlex.split(os.environ.get("MPIRUN", ""))


def definition(path, name, rows, count=None):
    text = [" ".join(map(str, row)) for row in rows]
    (path / name).write_text("====\nNItems {}\n====\n====\n====\n{}\n".format(
        len(rows) if count is None else count, "\n".join(text)))


def run(path, label, symmetry=True, layout="replicated", fail=None):
    command = (MPI if symmetry else []) + [HPHI, "-e", "sym.def" if symmetry else "raw.def"]
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    log = path / (label + ".log")
    with log.open("w") as stream:
        result = subprocess.run(command, cwd=str(path), env=env, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=120)
    text = log.read_text()
    if fail is not None:
        assert result.returncode != 0 and fail in text, text
        return None
    assert result.returncode == 0, "{}\n{}".format(command, text)
    if not symmetry and (path / "output/Eigenvalue.dat").exists():
        return np.sort(np.loadtxt(path / "output/Eigenvalue.dat")[:, 1])
    values = [float(line.split()[1]) for line in (path / "output/zvo_energy.dat").read_text().splitlines()
              if line.lstrip().startswith("Energy ")]
    return np.sort(values)


def tensor_operator(width, site, matrix, fermion=False):
    # Site 0 is the least significant tensor factor. Jordan-Wigner strings
    # precede it in the occupation ordering, including both Hubbard spins.
    result = np.ones((1, 1))
    for pos in reversed(range(width)):
        factor = matrix if pos == site else (np.diag([1, -1]) if fermion and pos < site else np.eye(2))
        result = np.kron(result, factor)
    return result


def prepare(model, length, nup=2, ndown=1, sector_test=None, sector_momenta=(0, 1)):
    path = ROOT / ("{}_up{}_down{}".format(model, nup, ndown) if model == "tJ" else model)
    if path.exists():
        shutil.rmtree(str(path))
    path.mkdir()
    width = length * (2 if model in ("Hubbard", "tJ") else 1)
    dim = 1 << width
    raw = np.zeros((dim, dim), dtype=complex)
    local = {}
    annihilator = {}
    if model != "Spin":
        for orbital in range(width):
            annihilator[orbital] = tensor_operator(width, orbital, np.array([[0, 1], [0, 0]]), True)

    def op(site, spin, other, other_spin):
        key = (site, spin, other, other_spin)
        if key not in local:
            if model == "Spin":
                assert site == other
                matrix = np.zeros((2, 2))
                matrix[spin, other_spin] = 1
                local[key] = tensor_operator(width, site, matrix)
            else:
                scale = 2 if model in ("Hubbard", "tJ") else 1
                local[key] = annihilator[scale*site+spin].T @ annihilator[scale*other+other_spin]
        return local[key]

    families = {}

    def add(family, row, matrix):
        nonlocal raw
        families.setdefault(family, []).append(row)
        raw += matrix

    def inter(indices, value):
        add("InterAll", list(indices) + [value.real, value.imag],
            value * (op(*indices[:4]) @ op(*indices[4:])))

    def pair(indices, value):
        inter(indices, value)
        # Dagger reverses both factor order and each creation/annihilation.
        a, s, b, t, c, u, d, v = indices
        inter((d, v, c, u, b, t, a, s), value.conjugate())

    for i in range(length):
        j, k = (i+1) % length, (i+2) % length
        spins = range(2 if model in ("Hubbard", "tJ") else 1)
        for s in spins:
            # Diagonal Transfer becomes EDChemi in the reader.
            add("Trans", [i, s, i, s, 0.17, 0], -0.17*op(i, s, i, s))
            if model != "Spin":
                z = 0.73 + 0.19j
                add("Trans", [i, s, j, s, z.real, z.imag], -z*op(i, s, j, s))
                add("Trans", [j, s, i, s, z.real, -z.imag], -z.conjugate()*op(j, s, i, s))
        if model == "Spin":
            exchange = op(i, 0, i, 1) @ op(j, 1, j, 0)
            add("Exchange", [i, j, 0.19], 0.19*(exchange + exchange.T))
            pair((i, 1, i, 0, j, 0, j, 1), 0.43+0.11j)
            inter((i, 0, i, 0, j, 0, j, 0), 0.31+0j)
            # Same-site product tests local matrix-unit contraction.
            inter((i, 1, i, 0, i, 0, i, 1), 0.23+0j)
        else:
            density_spin = 1 if model in ("Hubbard", "tJ") else 0
            pair((i, 0, j, 0, k, density_spin, k, density_spin), 0.21+0.07j)
            inter((i, 0, i, 0, j, density_spin, j, density_spin), 0.31+0j)
            # c_i^dag c_j c_j^dag c_i generates both density and -density-density.
            inter((i, 0, j, 0, j, 0, i, 0), 0.23+0j)
        if model in ("Hubbard", "tJ"):
            ni = op(i, 0, i, 0)+op(i, 1, i, 1)
            nj = op(j, 0, j, 0)+op(j, 1, j, 1)
            add("CoulombIntra", [i, 1.2], 1.2*(op(i, 0, i, 0) @ op(i, 1, i, 1)))
            add("CoulombInter", [i, j, 0.37], 0.37*(ni @ nj))
            add("Hund", [i, j, 0.29], -0.29*sum(op(i, s, i, s) @ op(j, s, j, s) for s in range(2)))
            exchange = op(i, 0, j, 0) @ op(j, 1, i, 1)
            add("Exchange", [i, j, 0.41], 0.41*(exchange + exchange.T))
            hopping = op(i, 0, j, 0) @ op(i, 1, j, 1)
            add("PairHop", [i, j, 0.13], 0.13*(hopping + hopping.T))

    if model in ("Hubbard", "tJ"):
        states = [s for s in range(dim) if sum((s >> (2*i)) & 1 for i in range(length)) == nup
                  and sum((s >> (2*i+1)) & 1 for i in range(length)) == ndown]
        if model == "tJ":
            states = [s for s in states if all(((s >> (2*i)) & 3) != 3 for i in range(length))]
        quantum = "Nup {}\nNdown {}\n".format(nup, ndown)
    else:
        states = [s for s in range(dim) if bin(s).count("1") == length//2]
        quantum = "2Sz 0\n" if model == "Spin" else "Ncond {}\n".format(length//2)
    raw = raw[np.ix_(states, states)]
    np.testing.assert_allclose(raw, raw.conj().T, atol=1e-12)
    definition(path, "loc.def", [[i, int(model == "Spin")] for i in range(length)],
               length if model == "Spin" else 0)
    namelist = "CalcMod calc.def\nModPara mod.def\nLocSpin loc.def\n"
    for family, rows in families.items():
        definition(path, family+".def", rows)
        namelist += "{} {}.def\n".format(family, family)
    (path / "raw.def").write_text(namelist)
    (path / "sym.def").write_text(namelist + "TransSym group.def\n")

    def settings(method, count):
        (path / "calc.def").write_text("CalcType {}\nCalcModel {}\nOutputMode 0\nOutputDataHead 1\n".format(
            method, {"Spin": 1, "SpinlessFermion": 7, "Hubbard": 0, "tJ": 9}[model]))
        (path / "mod.def").write_text("====\nModel_Parameters 0\n====\n====\n====\nCDataFileHead zvo\nCParaFileHead zqp\n====\n"
            "Nsite {}\n{}Lanczos_max 400\ninitial_iv -1\nexct {}\nLanczosEps 12\nLargeValue 100\nPreCG 0\n".format(length, quantum, count))

    settings(2, len(states))
    if model == "SpinlessFermion":
        # Raw spinless has no off-diagonal InterAll path. Its exactly equivalent
        # fully polarized Hubbard sector provides an independent raw reference.
        (path / "calc.def").write_text((path / "calc.def").read_text().replace("CalcModel 7", "CalcModel 0"))
        (path / "mod.def").write_text((path / "mod.def").read_text().replace(
            "Ncond {}".format(length//2), "Nup {}\nNdown 0".format(length//2)))
    if sector_test is None:
        np.testing.assert_allclose(run(path, "raw", symmetry=False), np.linalg.eigvalsh(raw), atol=2e-8, rtol=0)
    seen_dims = 0
    sector_spectra = []
    for momentum in range(length):
        if sector_test is not None and momentum not in sector_momenta:
            continue
        projector = np.zeros_like(raw)
        rows = []
        for g in range(length):
            ch = cmath.exp(-2j*np.pi*momentum*g/length)
            rows.append([g, ch.real, ch.imag])
            for col, state in enumerate(states):
                occupied = [i for i in range(width) if (state >> i) & 1]
                moved = [2*((i//2+g) % length)+i % 2 if model in ("Hubbard", "tJ") else (i+g) % length for i in occupied]
                sign = 1 if model == "Spin" else (-1)**sum(a > b for idx, a in enumerate(moved) for b in moved[idx+1:])
                target = sum(1 << i for i in moved)
                projector[states.index(target), col] += ch.conjugate()*sign/length
        for g in range(length):
            rows.extend([[g, i, (i+g) % length, 1] for i in range(length)])
        definition(path, "group.def", rows, length)
        text = (path / "group.def").read_text().replace("NItems", "NQPTrans")
        (path / "group.def").write_text(text)
        np.testing.assert_allclose(raw @ projector, projector @ raw, atol=1e-12)
        weights, vectors = np.linalg.eigh(projector)
        basis = vectors[:, weights > 0.5]
        dim_sector = basis.shape[1]
        seen_dims += dim_sector
        expected = np.linalg.eigvalsh(basis.conj().T @ raw @ basis)
        sector_spectra.append(expected)
        settings(3, dim_sector)
        if sector_test is not None:
            sector_test(path, model, length, momentum, states, raw, projector)
            continue
        for layout in ["replicated", "distributed"]:
            result = run(path, "k{}_{}".format(momentum, layout), layout=layout)
            np.testing.assert_allclose(result, expected, atol=3e-8, rtol=0)
            if model == "tJ":
                record = dict(line.split("=", 1) for line in (path / "output/zvo_symmetry_sector.dat").read_text().splitlines())
                assert record["model"] == "tJ" and int(record["full_dim"]) == len(states)
                assert int(record["sector_dim"]) == dim_sector
                assert int(record["fixed_nup"]) == nup and int(record["fixed_ndown"]) == ndown
                doublons = [float(line.split()[1]) for line in (path / "output/zvo_energy.dat").read_text().splitlines()
                            if line.lstrip().startswith("Doublon ")]
                assert len(doublons) == dim_sector and max(map(abs, doublons)) < 1e-12
    if sector_test is not None:
        return
    assert seen_dims == len(states)
    if model == "SpinlessFermion":
        assert np.max(np.abs(sector_spectra[1] - sector_spectra[-1])) > 0.1
    if model == "tJ":
        # These families vanish after projection, even with unequal couplings.
        for family in ["CoulombIntra", "PairHop"]:
            families[family][0][-1] += 0.6
            definition(path, family+".def", families[family])
        np.testing.assert_allclose(run(path, "projected_zero_terms"), expected, atol=3e-8, rtol=0)
    # Break translation invariance only in a newly supported diagonal family.
    rows = families["Trans"]
    rows[0][-2] += 0.5
    definition(path, "Trans.def", rows)
    run(path, "noninvariant", fail="Hamiltonian invariance failed")
    # Cancel the local perturbation through a different family. A validator
    # that compares families separately, or misses n_i^2 = n_i, rejects this.
    extra = [0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0]
    definition(path, "InterAll.def", families["InterAll"] + [extra])
    np.testing.assert_allclose(run(path, "cross_family_cancellation"), expected, atol=3e-8, rtol=0)
    if model == "Spin":
        # Invariant uniform transverse field must be rejected for fixed 2Sz.
        rows[0][-2] -= 0.5
        transverse = [[i, s, i, 1-s, 0.4, 0] for i in range(length) for s in range(2)]
        definition(path, "Trans.def", rows + transverse)
        definition(path, "InterAll.def", families["InterAll"])
        run(path, "not_sz_conserving", fail="does not conserve fixed Sz/Nup/Ndown")
    print("{}: {} raw levels and all {} momentum sectors verified".format(model, len(states), length))


if __name__ == "__main__":
    for model, length in [("Spin", 6), ("SpinlessFermion", 6), ("Hubbard", 4), ("tJ", 4)]:
        prepare(model, length)
    prepare("tJ", 4, 1, 1)
    prepare("tJ", 4, 2, 2)

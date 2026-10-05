"""Sector mTPQ against independent dense evolution, for every output column.

Only random initial vectors come from the test helper using HPhi's existing
RNG. Hamiltonians, signed projectors, propagation and observables are built
independently in Python. No statistical tolerances or thermal averaging.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import symmetry_general_terms as fixture

fixture.ROOT = Path("symmetry_tpq")
fixture.ROOT.mkdir(exist_ok=True)
INITIAL = str(Path(sys.argv[2]).resolve())


def invoke(path, command, label, env=None, failure=None):
    log = path / (label + ".log")
    with log.open("w") as stream:
        result = subprocess.run(command, cwd=str(path), env=env,
                                stdout=stream, stderr=subprocess.STDOUT, timeout=120)
    text = log.read_text()
    if failure:
        assert result.returncode != 0 and failure in text, text
    else:
        assert result.returncode == 0, text
    return text


def check_sector(path, model, length, momentum, states, raw, projector):
    # Choose P|r> for the smallest integer state in each orbit. Its coefficient
    # at r is real positive; this fixes the same mathematical gauge as HPhi.
    columns, covered = [], set()
    for index in range(len(states)):
        if index in covered or np.linalg.norm(projector[:, index]) < 1e-12:
            continue
        col = projector[:, index]
        covered.update(np.flatnonzero(abs(col) > 1e-12))
        columns.append(col / np.linalg.norm(col))
    basis = np.column_stack(columns)
    dim = basis.shape[1]
    np.testing.assert_allclose(basis.conj().T @ basis, np.eye(dim), atol=1e-12)
    h = basis.conj().T @ raw @ basis
    dtype = {"Spin": -1, "SpinlessFermion": 1, "Hubbard": 0, "tJ": -1}[model]
    calc = (path / "calc.def").read_text().replace("CalcType 3", "CalcType 1")
    calc += "InitialVecType {}\n".format(dtype)
    mod = (path / "mod.def").read_text().replace("Lanczos_max 400", "Lanczos_max 5")
    mod = mod.replace("LargeValue 100", "LargeValue 4").replace("initial_iv -1", "initial_iv 7")
    # exct is not a TPQ parameter and must not trigger the eigensolver guard.
    mod = mod.replace("exct {}".format(dim), "exct {}".format(dim + 1))
    (path / "mod.def").write_text(mod + "NumAve 2\n")
    for old in path.glob("initial_*_rank*.dat"):
        old.unlink()
    invoke(path, fixture.MPI + [INITIAL, "--dump-initial", str(dim), str(dtype), "7", "2"],
           "initial_k{}".format(momentum))
    vectors, first_norms = [], []
    ranks = sorted(path.glob("initial_0_rank*.dat"),
                   key=lambda p: int(p.stem.split("rank")[1]))
    for sample in range(2):
        values = []
        for rank in range(len(ranks)):
            lines = (path / "initial_{}_rank{}.dat".format(sample, rank)).read_text().splitlines()
            first = float(lines[0])
            values.extend(complex(*map(float, line.split())) for line in lines[1:])
        assert len(values) == dim
        vectors.append(np.array(values))
        first_norms.append(first)
    # Local diagonal observables are invariant over each group orbit.
    number = length if model == "Spin" else (length // 2 if model == "SpinlessFermion" else 3)
    sz = 0.5 if model in ("Hubbard", "tJ") else 0.0
    doublon_raw = np.array([sum(((s >> (2*i)) & 3) == 3 for i in range(length))
                           if model in ("Hubbard", "tJ") else 0 for s in states])
    doublon = basis.conj().T @ (doublon_raw[:, None] * basis)
    doublon2 = basis.conj().T @ ((doublon_raw**2)[:, None] * basis)
    expected = []
    for sample, vector in enumerate(vectors):
        ss, norms, flct = [], [], []
        norm = first_norms[sample]
        np.testing.assert_allclose(np.vdot(vector, vector), 1, atol=1e-12)
        for step in range(5):
            hv = h @ vector
            energy = np.vdot(vector, hv).real
            energy2 = np.vdot(hv, hv).real
            d = np.vdot(vector, doublon @ vector).real
            d2 = np.vdot(vector, doublon2 @ vector).real
            beta = 2*step/(4*length-energy)
            ss.append([beta, energy, energy2, d, number, step])
            norms.append([beta, norm, first_norms[sample], step])
            if step:
                flct.append([beta, number, number**2, d, d2, sz, sz**2, step])
            vector = 4*vector - hv/length
            norm = np.linalg.norm(vector)
            vector /= norm
        expected.append([np.array(ss), np.array(norms), np.array(flct)])
    for layout in ["default", "replicated"]:
        env = dict(os.environ)
        env.pop("HPHI_SYMMETRY_BASIS_LAYOUT", None)
        if layout != "default":
            env["HPHI_SYMMETRY_BASIS_LAYOUT"] = layout
        for aggregate in [0, 1]:
            (path / "calc.def").write_text(calc + "OutputGreenFormat {}\n".format(aggregate))
            if (path / "output").exists():
                shutil.rmtree(path / "output")
            label = "k{}_{}_aggregate{}".format(momentum, layout, aggregate)
            text = invoke(path, fixture.MPI + [fixture.HPHI, "-e", "sym.def"], label, env)
            if layout == "default":
                assert "distributed (default for TransSym TPQ)" in text
            record = dict(line.split("=", 1) for line in (path / "output/zvo_symmetry_sector.dat").read_text().splitlines())
            assert record["ensemble"] == "single_symmetry_sector"
            assert record["calc_type"] == "TPQ" and int(record["sector_dim"]) == dim
            assert int(record["num_ave"]) == 2 and float(record["large_value"]) == 4
            assert int(record["initial_vec_type"]) == dtype and int(record["mpi_ranks"]) == len(ranks)
            for sample in range(2):
                for family, reference in zip(["SS", "Norm", "Flct"], expected[sample]):
                    name = "{}_tpq.dat".format(family) if aggregate else "{}_rand{}.dat".format(family, sample)
                    actual = np.loadtxt(path / "output" / ("zvo_" + name), ndmin=2)
                    if aggregate:
                        actual = actual[actual[:, 0] == sample]
                        actual = np.column_stack((actual[:, 2:], actual[:, 1]))
                    np.testing.assert_allclose(actual, reference, rtol=2e-11, atol=2e-11,
                                               err_msg="{} {} {} sample{}".format(model, label, family, sample))
            # Preserve outputs for inspection instead of only the last run.
            shutil.copytree(path / "output", path / label)
    # Feature gates must report TPQ-specific scope, before any MPI matvec.
    if model == "Spin" and momentum == 0:
        (path / "calc.def").write_text(calc)
        for option in ["ReStart 1", "InputEigenVec 1", "OutputEigenVec 1", "InputHam 1", "OutputHam 1"]:
            (path / "calc.def").write_text(calc + option + "\n")
            invoke(path, [fixture.HPHI, "-e", "sym.def"], option.split()[0],
                   failure=("OutputHam is only defined for FullDiag" if option.startswith("OutputHam")
                            else "does not support " + option.split()[0]))
        (path / "calc.def").write_text(calc)
    print("{} k={}: {} dimensions, all TPQ columns match independent evolution".format(model, momentum, dim))


for model, length in [("Spin", 6), ("SpinlessFermion", 6), ("Hubbard", 4), ("tJ", 4)]:
    fixture.prepare(model, length, sector_test=check_sector)

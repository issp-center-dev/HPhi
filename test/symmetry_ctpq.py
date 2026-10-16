"""Sector cTPQ against independent dense evolution, for every output column.

Only random initial vectors come from the test helper using HPhi's existing
RNG. Hamiltonians, signed projectors, propagation and observables are built
independently in Python. No statistical tolerances or thermal averaging.
"""
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import symmetry_correlation as correlation
import symmetry_general_terms as fixture

fixture.ROOT = Path("symmetry_ctpq")
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


def check_schedule(path, model, length, momentum, states, raw, projector, explicit):
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
    requests = correlation.step_requests(model, length)
    operators = correlation.sector_operators(model, length, states, basis, requests)
    dtype = {"Spin": -1, "SpinlessFermion": 1, "Hubbard": 0, "tJ": -1}[model]
    calc = (path / "calc.def").read_text().replace("CalcType 3", "CalcType 5")
    calc += "InitialVecType {}\n".format(dtype)
    mod = (path / "mod.def").read_text().replace("Lanczos_max 400", "Lanczos_max 5")
    mod = mod.replace("LargeValue 100", "LargeValue 4").replace("initial_iv -1", "initial_iv 7")
    # exct is not a TPQ parameter and must not trigger the eigensolver guard.
    mod = mod.replace("exct {}".format(dim), "exct {}".format(dim + 1))
    # The general fixture is reused for two schedules. Remove our previous
    # additions before writing each one.
    mod = "\n".join(line for line in mod.splitlines() if not line.startswith(("NumAve", "ExpandCoef"))) + "\n"
    base_names = "\n".join(line for line in (path / "sym.def").read_text().splitlines()
                           if not line.startswith("InvTemp")) + "\n"
    order = 10 if model in ("Spin", "tJ") else 4
    betas = np.array([0, .13, .29, .29, .52]) if explicit else np.arange(5)*.25
    orders = [2, 7, 3, 6, 5] if explicit else [order]*5
    (path / "sym.def").write_text(base_names + ("InvTemp beta.def\n" if explicit else ""))
    correlation.add_requests(path, "sym.def", requests)
    (path / "beta.def").write_text("".join("{} {} {} 0\n".format(b, n, i % 2)
                                            for i, (b, n) in enumerate(zip(betas, orders))))
    (path / "mod.def").write_text(mod + "NumAve 2\n" +
                                ("ExpandCoef 4\n" if order == 4 and not explicit else ""))
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
    expected, vectors_by_step = [], []
    for sample, vector in enumerate(vectors):
        ss, norms, flct, step_vectors = [], [], [], []
        norm = first_norms[sample]
        np.testing.assert_allclose(np.vdot(vector, vector), 1, atol=1e-12)
        for step in range(5):
            step_vectors.append(vector.copy())
            hv = h @ vector
            energy = np.vdot(vector, hv).real
            energy2 = np.vdot(hv, hv).real
            d = np.vdot(vector, doublon @ vector).real
            d2 = np.vdot(vector, doublon2 @ vector).real
            beta = betas[step]
            ss.append([beta, energy, energy2, d, number, step])
            norms.append([beta, norm, first_norms[sample], step])
            flct.append([beta, number, number**2, d, d2, sz, sz**2, step])
            if step < 4:
                # Evaluate the polynomial independently in the exact eigenbasis,
                # instead of mirroring HPhi's repeated sparse matrix products.
                eig, rotation = np.linalg.eigh(h)
                x = -.5*(betas[step+1]-betas[step])*eig
                polynomial = sum(x**n / math.factorial(n)
                                 for n in range(orders[step]+1))
                vector = rotation @ (polynomial * (rotation.conj().T @ vector))
                norm = np.linalg.norm(vector)
                vector /= norm
        expected.append([np.array(ss), np.array(norms), np.array(flct)])
        vectors_by_step.append(step_vectors)
    correlation.assert_nonzero_operators(
        operators, [vector for sample in vectors_by_step for vector in sample]
    )
    emitted_steps = [0, 1, 3] if explicit else list(range(5))
    for layout in ["default", "replicated"]:
        env = dict(os.environ)
        env.pop("HPHI_SYMMETRY_BASIS_LAYOUT", None)
        if layout != "default":
            env["HPHI_SYMMETRY_BASIS_LAYOUT"] = layout
        for aggregate in [0, 1]:
            (path / "calc.def").write_text(calc + "OutputGreenFormat {}\n".format(aggregate))
            if (path / "output").exists():
                shutil.rmtree(path / "output")
            label = "k{}_schedule{}_{}_aggregate{}".format(momentum, int(explicit), layout, aggregate)
            text = invoke(path, fixture.MPI + [fixture.HPHI, "-e", "sym.def"], label, env)
            if layout == "default":
                assert "distributed (default for TransSym cTPQ)" in text
            record = dict(line.split("=", 1) for line in (path / "output/zvo_symmetry_sector.dat").read_text().splitlines())
            assert record["ensemble"] == "single_symmetry_sector"
            assert record["calc_type"] == "cTPQ" and int(record["sector_dim"]) == dim
            assert int(record["num_ave"]) == 2 and float(record["large_value"]) == 4
            assert int(record["initial_vec_type"]) == dtype and int(record["mpi_ranks"]) == len(ranks)
            assert int(record["canonical_tpq_steps"]) == 5
            assert record["beta_schedule"] == ("explicit" if explicit else "uniform")
            if explicit:
                for row in range(5):
                    np.testing.assert_allclose(list(map(float, record["invtemp_row_{}".format(row)].split())),
                                               [betas[row], orders[row]], atol=1e-15)
            else:
                assert float(record["beta_step"]) == .25 and int(record["expand_coef"]) == order
            for sample in range(2):
                for family, reference in zip(["SS", "Norm", "Flct"], expected[sample]):
                    name = "{}_tpq.dat".format(family) if aggregate else "{}_rand{}.dat".format(family, sample)
                    actual = np.loadtxt(path / "output" / ("zvo_" + name), ndmin=2)
                    if aggregate:
                        actual = actual[actual[:, 0] == sample]
                        actual = np.column_stack((actual[:, 2:], actual[:, 1]))
                    np.testing.assert_allclose(actual, reference, rtol=2e-11, atol=2e-11,
                                               err_msg="{} {} {} sample{}".format(model, label, family, sample))
                for step in range(5):
                    if step not in emitted_steps:
                        if not aggregate:
                            for kind in requests:
                                missing = path / "output" / "zvo_{}_set{}step{}.dat".format(
                                    correlation.FILES[kind], sample, step
                                )
                                assert not missing.exists(), (label, missing)
                        continue
                    if aggregate:
                        name = "zvo_{}_tpq.dat"
                        flt = lambda data, s=sample, n=step: data[
                            (data[:, 0] == s) & (data[:, 1] == n)
                        ]
                    else:
                        name = "zvo_{{}}_set{}step{}.dat".format(sample, step)
                        flt = None
                    correlation.check_step_file(
                        path, name, requests, operators,
                        vectors_by_step[sample][step],
                        prefix=2 if aggregate else 0,
                        prefix_filter=flt,
                        label="{} {} sample{} step{}".format(model, label, sample, step),
                    )
            if aggregate:
                correlation.check_aggregate_index(
                    path, "zvo_{}_tpq.dat", requests,
                    [(sample, step) for sample in range(2)
                     for step in emitted_steps],
                    label=label,
                )
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
        if explicit:
            original_beta = (path / "beta.def").read_text()
            for row in (0, 2, 4):
                lines = original_beta.splitlines()
                fields = lines[row].split()
                fields[-1] = "1"
                lines[row] = " ".join(fields)
                (path / "beta.def").write_text("\n".join(lines) + "\n")
                (path / "calc.def").write_text(calc)
                invoke(path, fixture.launcher(True) + [fixture.HPHI, "-e", "sym.def"], "reject_eigen_{}".format(row),
                       failure="does not support InvTemp eigenvector output")
                assert not list((path / "output").glob("*eigen*"))
            (path / "beta.def").write_text(original_beta)
    print("{} k={}: {} dimensions, cTPQ schedule {} matches independent evolution".format(model, momentum, dim, explicit))


def check_sector(path, model, length, momentum, states, raw, projector):
    original_mod = (path / "mod.def").read_text()
    original_calc = (path / "calc.def").read_text()
    for explicit in (False, True):
        (path / "mod.def").write_text(original_mod)
        (path / "calc.def").write_text(original_calc)
        check_schedule(path, model, length, momentum, states, raw, projector, explicit)


for model, length in [("Spin", 6), ("SpinlessFermion", 6), ("Hubbard", 4), ("tJ", 4)]:
    fixture.prepare(model, length, sector_test=check_sector)

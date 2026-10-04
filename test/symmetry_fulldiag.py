"""All sector eigenvalues: independent tensor H, raw FullDiag, LAPACK and MPI backends."""
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import symmetry_general_terms as fixture

fixture.ROOT = Path("symmetry_fulldiag")
fixture.ROOT.mkdir(exist_ok=True)
spectra = []


def execute(path, label, calc, launcher=(), layout=None, failure=None):
    (path / "calc.def").write_text(calc)
    env = dict(os.environ)
    env.pop("HPHI_SYMMETRY_BASIS_LAYOUT", None)
    if layout is not None:
        env["HPHI_SYMMETRY_BASIS_LAYOUT"] = layout
    log = path / (label + ".log")
    with log.open("w") as stream:
        result = subprocess.run(list(launcher) + [fixture.HPHI, "-e", "sym.def"],
                                cwd=str(path), env=env, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=60)
    text = log.read_text()
    if failure:
        assert result.returncode != 0 and failure in text, text
        return text
    assert result.returncode == 0, text
    return text


def check_sector(path, model, length, momentum, states, raw, projector):
    weights, vectors = np.linalg.eigh(projector)
    basis = vectors[:, weights > 0.5]
    expected = np.linalg.eigvalsh(basis.conj().T @ raw @ basis)
    dim = len(expected)
    calc = (path / "calc.def").read_text().replace("CalcType 3", "CalcType 2")
    mod = (path / "mod.def").read_text()
    if momentum == 0:
        spectra.clear()
        # Use raw canonical FullDiag only in serial: its site decomposition
        # is unrelated to the sector solver's column ownership.
        rawcalc = calc
        if model == "SpinlessFermion":
            rawcalc = calc.replace("CalcModel 7", "CalcModel 0")
            (path / "mod.def").write_text(mod.replace("Ncond {}".format(length//2),
                                                       "Nup {}\nNdown 0".format(length//2)))
        (path / "calc.def").write_text(rawcalc + "Solver 0\n")
        np.testing.assert_allclose(fixture.run(path, "raw_fulldiag", symmetry=False),
                                   np.linalg.eigvalsh(raw), rtol=0, atol=2e-8)
        (path / "mod.def").write_text(mod)
    backends = [(0, [])]
    if fixture.MPI and os.environ.get("HPHI_HAS_SCALAPACK") == "1":
        backends.append((1, fixture.MPI))
    if os.environ.get("HPHI_HAS_ELPA") == "1":
        backends.append((3, fixture.MPI))
    for solver, launcher in backends:
        if (path / "output").exists():
            shutil.rmtree(path / "output")
        settings = calc + "Solver {}\nNGPU 0\n".format(solver)
        label = "k{}_solver{}".format(momentum, solver)
        # ELPA explicitly rejects a sector smaller than its process grid.
        ranks = 1
        for flag in ("-np", "-n", "--np", "--n"):
            if flag in launcher:
                ranks = int(launcher[launcher.index(flag)+1])
                break
        divisors = [d for d in range(1, int(ranks**0.5)+1) if ranks % d == 0]
        grid_max = ranks // max(divisors)
        if solver == 3 and dim < grid_max:
            execute(path, label, settings, launcher, failure="smaller than the process grid")
            assert not (path / "output/zvo_energy_sector.dat").exists()
            continue
        text = execute(path, label, settings, launcher)
        assert "replicated (default for TransSym FullDiag)" in text
        actual = np.loadtxt(path / "output/zvo_energy_sector.dat", ndmin=2)
        np.testing.assert_array_equal(actual[:, 0], np.arange(dim))
        np.testing.assert_allclose(actual[:, 1], expected, rtol=0, atol=3e-8,
                                   err_msg="{} {}".format(model, label))
        assert not (path / "output/Eigenvalue.dat").exists()
        assert not list((path / "output").glob("*phys*"))
        assert "Symmetry matvec plan:" not in text
        assert "Symmetry FullDiag storage: dimension={} metadata=replicated".format(dim) in text
        assert "raw_basis_elements=0" in text
        if solver == 1 and ranks > 1:
            assert "matrix=replicated eigenvectors=distributed" in text
        if solver == 3 and ranks > 1:
            assert "matrix=column_panel eigenvectors=distributed" in text
        manifest = dict(line.split("=", 1) for line in (path / "output/zvo_symmetry_sector.dat").read_text().splitlines())
        assert manifest["calc_type"] == "FullDiag" and int(manifest["sector_dim"]) == dim
        assert int(manifest["full_dim"]) == len(states) and manifest["basis_layout"] == "replicated"
        assert manifest["output_scope"] == "eigenvalues" and int(manifest["solver_id"]) == solver
        assert manifest["eigenvalue_file"] == "zvo_energy_sector.dat"
        if solver == 0:
            spectra.extend(actual[:, 1])
        shutil.copytree(path / "output", path / label)
    if momentum == length-1:
        np.testing.assert_allclose(sorted(spectra), np.linalg.eigvalsh(raw), rtol=0, atol=3e-8)
    if model == "Spin" and momentum == 0:
        execute(path, "layout_reject", calc, layout="distributed",
                failure="distributed symmetry basis is supported for TransSym Lanczos, TPQ, CG, TimeEvolution and cTPQ runs only")
        for option in ["ReStart", "InputEigenVec", "OutputEigenVec", "InputHam", "OutputHam"]:
            execute(path, "reject_"+option, calc + option + " 1\n",
                    failure="does not support " + option)
        fixture.definition(path, "one.def", [[0, 0, 0, 0]])
        original = (path / "sym.def").read_text()
        (path / "sym.def").write_text(original + "OneBodyG one.def\n")
        execute(path, "correlation_reject", calc, failure="sector FullDiag outputs eigenvalues only")
        (path / "sym.def").write_text(original)
        for solver, launcher in backends:
            if solver == 3 and dim < grid_max:
                continue
            output = path / "output/zvo_energy_sector.dat"
            if output.exists():
                output.unlink()
            output.mkdir()
            # Verify any output-open failure returns a nonzero code without
            # assuming the platform-specific strerror text.
            (path / "calc.def").write_text(calc + "Solver {}\nNGPU 0\n".format(solver))
            with (path / "write_failure_{}.log".format(solver)).open("w") as stream:
                result = subprocess.run(launcher + [fixture.HPHI, "-e", "sym.def"], cwd=str(path),
                                        env=dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT="replicated"),
                                        stdout=stream, stderr=subprocess.STDOUT, timeout=30)
            assert result.returncode != 0
            output.rmdir()
    print("{} k={}: {} sector eigenvalues match; tested solvers {}".format(model, momentum, dim,
                                                                             [s for s, _ in backends]))


for model, length in [("Spin", 6), ("SpinlessFermion", 6), ("Hubbard", 4), ("tJ", 4)]:
    fixture.prepare(model, length, sector_test=check_sector, sector_momenta=range(length))

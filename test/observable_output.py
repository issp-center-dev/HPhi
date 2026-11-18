"""Regressions for reported tJ Sz and spinless particle number.

Lanczos is used only to exercise its Sz output path. Independent Fock-space
matrices check the energy and saved-vector residual; other cases use CG.
"""
from pathlib import Path
import os
import re
import shlex
import struct
import subprocess
import sys
import tempfile

import numpy as np
from symmetry_spingc_common import definition, join_rank_vectors, mpi_size
from symmetry_kondo_common import read_vector_parts
from symmetry_nconserved_common import expert_case as tj_input
from symmetry_nconserved_reference import (
    Case as TJCase, DEFAULT, raw_words as tj_words, sector_basis as tj_basis,
    operator_matrix, physical_terms)
from symmetry_fermion_gc_common import expert_case as gc_input, join_gc_vectors
from symmetry_fermion_gc_reference import (
    Case as GCCase, make_reference, pair_rows)

HPHI = str(Path(sys.argv[1]).resolve())
MODE = sys.argv[2]
PROBE = str(Path(sys.argv[3]).resolve()) if len(sys.argv) > 3 else HPHI
MPI = shlex.split(os.environ.get("MPIRUN", ""))
LAYOUT = os.environ.get("HPHI_SYMMETRY_BASIS_LAYOUT", "distributed")
Path(f"{MODE}_output").mkdir(exist_ok=True)
# Never accept vector/output files left behind by a previous invocation.
ROOT = Path(tempfile.mkdtemp(prefix=f"np{mpi_size(MPI)}_{LAYOUT}_",
                             dir=f"{MODE}_output"))
ENV = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=LAYOUT)


def run(path, sector, capture=False):
    env = ENV.copy()
    if not sector:
        env.pop("HPHI_SYMMETRY_BASIS_LAYOUT", None)
    if capture:
        env.update(HPHI_TEST_SYMMETRY_CAPTURE="1", HPHI_TEST_SYMMETRY_ACTION="lanczos")
    with (path / "run.log").open("w") as log:
        result = subprocess.run(MPI + [PROBE if capture else HPHI, "-e", "sym.def"], cwd=path,
                                env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=120)
    text = (path / "run.log").read_text(errors="replace")
    assert result.returncode == 0, f"{path.name}: rc={result.returncode}\n{text}"
    return text


def raw_input(path):
    p = path / "sym.def"
    p.write_text("".join(line for line in p.read_text().splitlines(True)
                         if not line.startswith("TransSym ")))


def files(path):
    return sorted((path / "output").glob("zvo_eigenvec_0_rank_*.dat"),
                  key=lambda p: int(p.stem.rsplit("_", 1)[1]))


def legacy_vector(path):
    parts = []
    for p in files(path):
        data = p.read_bytes()
        step, dim = struct.unpack("=iQ", data[:12])
        assert step >= 0 and len(data) == 12 + 16 * (dim + 1), p
        parts.append(np.frombuffer(data[12:], dtype=np.complex128)[1:])
    assert parts, path
    return np.concatenate(parts)


def check_state(path, vector, h, tolerance):
    assert vector.shape == (len(h),) and np.isfinite(vector).all(), path
    np.testing.assert_allclose(np.vdot(vector, vector), 1., atol=1e-10, rtol=0)
    eigenvalues = np.linalg.eigvalsh(h)
    energy = float(np.vdot(vector, h @ vector).real)
    scale = max(1., float(np.max(np.abs(eigenvalues))))
    residual = float(np.linalg.norm(h @ vector - energy * vector) / scale)
    assert abs(energy - eigenvalues[0]) < tolerance * scale, (path, energy, eigenvalues[0])
    assert residual < tolerance, (path, residual)
    energy_text = (path / "output/zvo_energy.dat").read_text()
    reported = float(re.search(r"^\s*Energy\s+(\S+)", energy_text, re.M)[1])
    assert abs(reported - energy) < 1e-9 * scale, (path, reported, energy)
    print(f"{path.name}: E={energy:.12g}, residual={residual:.3e}", flush=True)
    return energy_text


def tj_case(twosz, sector):
    # Eight sites also allow the raw path in the baseline np16 job.
    case = TJCase("tJNConserved", 8, 2, 1)
    parameters = dict(DEFAULT, soc=0., soc2=0., hx=0., w=0.)
    # Restrict the independent N-only oracle to the canonical 2Sz block.
    words = tj_words(case)
    mask = np.array([sum(1 if o % 2 == 0 else -1
                         for o in range(2 * case.sites) if int(w) & (1 << o)) == twosz
                     for w in words])
    words = words[mask]
    h = operator_matrix(case, words, physical_terms(case.sites, parameters))
    if sector:
        basis, _ = tj_basis(case, words)
        h = basis.conj().T @ h @ basis
    path = tj_input(ROOT / f"sz{twosz}_{'sector' if sector else 'raw'}", case,
                    method=0, parameters=parameters,
                    options={"CalcMod": {"CalcEigenVec": 1, "OutputEigenVec": int(not sector)},
                             "ModPara": {"2Sz": twosz, "LanczosEps": 14}})
    if not sector:
        raw_input(path)
    log = run(path, sector, capture=sector)
    # Sector Lanczos forbids public vector I/O; use the existing test-only
    # capture hook, which shares the production sources and compiler flags.
    vector = (read_vector_parts(path, "sector_final_sample0_step0.rank", len(h))
              if sector else legacy_vector(path))
    text = check_state(path, vector, h, 1e-7)
    assert "Accuracy check" in log, log
    actual = float(re.search(r"^\s*Sz\s+(\S+)", text, re.M)[1])
    print(f"{path.name}: reported Sz={actual}, expected={twosz / 2}", flush=True)
    assert abs(actual - twosz / 2) < 1e-12, (path, actual, twosz)


def spinless_case(name, sector, *, canonical=False, parity=0, pairs=False):
    case = GCCase("SpinlessFermionGC", 8, parity, 1)
    pairing = pair_rows(case.model, case.sites) if pairs else None
    ref = make_reference(case, pairs=pairing)
    words, h, basis = ref.words, ref.hamiltonian, ref.basis
    options = {"CalcMod": {"OutputEigenVec": 1}, "ModPara": {"LanczosEps": 22}}
    if canonical:
        options["CalcMod"]["CalcModel"] = 7
        options["ModPara"]["Ncond"] = 3
        mask = np.array([int(w).bit_count() == 3 for w in words])
        words = words[mask]
        h = h[np.ix_(mask, mask)]
        restricted = basis[mask, :]
        basis = restricted[:, np.linalg.norm(restricted, axis=0) > .5]
    if not sector:
        basis = np.eye(len(words))
    matrix = basis.conj().T @ h @ basis
    path = gc_input(ROOT / name, case, method=3, options=options, pairs=pairing)
    if not sector:
        raw_input(path)
    definition(path, "density.def", [(j, 0, j, 0) for j in range(case.sites)])
    with (path / "sym.def").open("a") as f:
        f.write("OneBodyG density.def\n")
    log = run(path, sector)
    if sector:
        reader = join_rank_vectors if canonical else join_gc_vectors
        vector = reader(files(path), len(matrix))
    else:
        vector = legacy_vector(path)
    check_state(path, vector, matrix, 1e-10)
    weights = abs(basis @ vector)**2
    counts = np.array([int(w).bit_count() for w in words])
    expected = float(weights @ counts)
    assert expected > .5, (name, expected)
    if pairs:
        assert abs(expected - round(expected)) > 1e-4, (name, expected)
    densities = np.loadtxt(path / "output/zvo_cisajs_eigen0.dat", ndmin=2)
    assert densities.shape == (case.sites, 6), densities
    saved = float(np.sum(densities[:, 4]))
    np.testing.assert_allclose(densities[:, 5], 0., atol=1e-10, rtol=0)
    assert abs(saved - expected) < 1e-8, (name, saved, expected)
    row = re.search(r"^\s*i=\s*0\s+Energy=\s*(\S+)\s+N=\s*(\S+)\s+Sz=\s*(\S+)", log, re.M)
    assert row, log
    actual, sz = float(row[2]), float(row[3])
    print(f"{name}: stdout N={actual:.6f}, vector N={expected:.12g}, "
          f"density sum={saved:.12g}", flush=True)
    assert abs(sz) < 1e-12, (name, sz)
    # phys() prints six digits after the decimal point.
    assert abs(actual - expected) <= 5.1e-7, (name, actual, expected)


if MODE == "tj_sz":
    cases = [(f"{s}/{sector}", lambda s=s, sector=sector: tj_case(s, sector))
             for s in (2, -2, 0) for sector in (False, True)]
elif MODE == "spinless_number":
    cases = [(name, lambda name=name, sector=sector, options=options:
              spinless_case(name, sector, **options))
             for name, sector, options in (
                 ("canonical_raw", False, {"canonical": True}),
                 ("canonical_sector", True, {"canonical": True}),
                 ("gc_raw", False, {}),
                 ("gc_sector", True, {}),
                 ("paired_even_sector", True, {"parity": 1, "pairs": True}),
                 ("paired_odd_sector", True, {"parity": -1, "pairs": True}))]
else:
    raise ValueError(MODE)
failures = []
for label, check in cases:
    try:
        check()
    except (AssertionError, subprocess.TimeoutExpired) as error:
        failures.append(label)
        print(f"FAIL {label}: {error}", flush=True)
assert not failures, failures
print(f"{MODE}: {len(cases)} cases passed ({LAYOUT})", flush=True)

"""Literal fermion strings are the oracle for InterAll contraction and TE.

FullDiag checks run serially with LAPACK. Identity-sector CG/TE checks use
MPIRUN and distributed vectors; a four-site raw GC case also tests site MPI.
"""
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile

import numpy as np


HPHI = str(Path(sys.argv[1]).resolve())
MPI = shlex.split(os.environ.get("MPIRUN", ""))
ROOT = Path(tempfile.mkdtemp(prefix="interall_contraction_", dir=".")).resolve()
ENV = dict(os.environ, OMP_NUM_THREADS="1",
           HPHI_SYMMETRY_BASIS_LAYOUT="distributed")


def definition(path, rows, count=None, label="NItems"):
    path.write_text("====\n{} {}\n====\n====\n====\n".format(
        label, len(rows) if count is None else count)
        + "".join(" ".join(map(str, row)) + "\n" for row in rows))


def pair(indices, value):
    sites = list(zip(indices[::2], indices[1::2]))
    dagger = tuple(x for site in reversed(sites) for x in site)
    return [(*indices, value.real, value.imag),
            (*dagger, value.real, -value.imag)]


def matrix(states, transfer, interall):
    """Apply the input strings right to left, before projecting the final state.

    This deliberately does not use the parser's reordered four-fermion form.
    """
    lookup = {word: i for i, word in enumerate(states)}
    h = np.zeros((len(states), len(states)), complex)
    for rows, prefactor in ((transfer, -1), (interall, 1)):
        for row in rows:
            indices = row[:-2]
            operators = [(i % 2 == 0, 2 * indices[2*i] + indices[2*i+1])
                         for i in range(len(indices)//2)]
            for col, word in enumerate(states):
                value = prefactor * complex(*row[-2:])
                for create, orbital in reversed(operators):
                    occupied = bool(word & (1 << orbital))
                    if occupied == create:
                        break
                    value *= (-1)**((word & ((1 << orbital)-1)).bit_count())
                    word ^= 1 << orbital
                else:
                    if word in lookup:
                        h[lookup[word], col] += value
    np.testing.assert_allclose(h, h.conj().T, atol=1e-13, rtol=0)
    return h


def basis(model, polarized=False, nsite=None):
    kondo = model.startswith("Kondo")
    length = nsite or (4 if kondo else 3)
    gc = model.endswith("GC")
    fixed_n = model.endswith("NConserved")
    # Four itinerant electrons avoid a separate low-filling underflow in the
    # legacy raw Kondo basis enumerator; the three-site contractions stay active.
    number = (4 if kondo else 1 if polarized else 2) + int(kondo)
    sz = 1 if polarized or kondo else 0
    words = []
    for word in range(1 << (2*length)):
        if kondo and (word & 3).bit_count() != 1:
            continue
        if model.startswith("tJ") and any((word >> (2*i)) & 3 == 3
                                          for i in range(length)):
            continue
        if not gc and word.bit_count() != number:
            continue
        if not gc and not fixed_n:
            actual_sz = sum(((word >> (2*i)) & 1) - ((word >> (2*i+1)) & 1)
                            for i in range(length))
            if actual_sz != sz:
                continue
        words.append(word)
    quantum = "" if gc else f"Ncond {number-int(kondo)}\n"
    if not gc and not fixed_n:
        quantum += f"2Sz {sz}\n"
    code = (5 if gc else 2) if kondo else (
        (10 if gc else 9) if model.startswith("tJ") else (3 if gc else 0))
    return length, code, quantum, words


def prepare(label, model, transfer, interall, *, sector, polarized=False, nsite=None):
    path = ROOT / label
    path.mkdir()
    length, code, quantum, words = basis(model, polarized, nsite)
    calc = f"CalcModel {code}\nOutputMode 0\nOutputDataHead 1\n"
    mod = ("====\nModel_Parameters 0\n====\n====\n====\n"
           "CDataFileHead zvo\nCParaFileHead zqp\n====\n"
           f"Nsite {length}\n{quantum}initial_iv -1\nexct 1\n"
           "LanczosEps 22\nLanczosTarget 1\nLargeValue 100\nPreCG 0\n")
    (path / "calc.def").write_text(calc + "CalcType 2\n")
    (path / "mod.def").write_text(mod + "Lanczos_max 400\n")
    definition(path / "loc.def",
               [(i, int(model.startswith("Kondo") and i == 0)) for i in range(length)],
               int(model.startswith("Kondo")))
    names = "CalcMod calc.def\nModPara mod.def\nLocSpin loc.def\n"
    for keyword, rows in (("Trans", transfer), ("InterAll", interall)):
        if rows:
            definition(path / (keyword + ".def"), rows)
            names += f"{keyword} {keyword}.def\n"
    if sector:
        definition(path / "group.def", [(0, 1, 0)]
                   + [(0, i, i, 1) for i in range(length)], 1, "NQPTrans")
        names += "TransSym group.def\n"
    (path / "namelist.def").write_text(names)
    return path, words, calc, mod, names


def run(path, label, sector, raw_mpi=False):
    command = (MPI if (sector or raw_mpi) and label != "fulldiag" else []) + [
        HPHI, "-e", "namelist.def"]
    env = dict(ENV, HPHI_SYMMETRY_BASIS_LAYOUT=(
        "distributed" if sector and label != "fulldiag" else "replicated"))
    with (path / (label + ".log")).open("w") as log:
        result = subprocess.run(command, cwd=path, env=env, stdout=log,
                                stderr=subprocess.STDOUT, timeout=300)
    assert result.returncode == 0, (path, (path / (label + ".log")).read_text())


def vector(path, stem, sector):
    pieces = []
    for file in sorted((path / "output").glob(stem + "_rank_*.dat"),
                       key=lambda p: int(p.stem.split("_rank_")[1])):
        data = file.read_bytes()
        if sector:
            version = struct.unpack_from("<Q", data, 8)[0]
            assert version in (1, 2), (file, version)
            header_size = 224 if version == 1 else 240
            pieces.extend(np.frombuffer(data[header_size:], dtype="<c16"))
        else:
            _, size = struct.unpack("=iQ", data[:12])
            part = np.frombuffer(data[12:], dtype=np.complex128)
            assert len(part) == size + 1
            pieces.extend(part[1:])
    assert pieces, (path, stem)
    return np.array(pieces)


def spectra(label, model, transfer, interall, *, sector, polarized=False, cg=False,
            nsite=None, raw_mpi=False):
    if model.startswith("Kondo"):
        # Keep the localized spin first, followed by three itinerant sites.
        def shift(rows):
            return [tuple(x+1 if i % 2 == 0 else x for i, x in enumerate(row[:-2]))
                    + tuple(row[-2:]) for row in rows]
        transfer, interall = shift(transfer), shift(interall)
    path, words, calc, mod, _ = prepare(
        label, model, transfer, interall, sector=sector, polarized=polarized, nsite=nsite)
    h = matrix(words, transfer, interall)
    run(path, "fulldiag", sector)
    output = "zvo_energy_sector.dat" if sector else "Eigenvalue.dat"
    actual = np.loadtxt(path / "output" / output, ndmin=2)[:, 1]
    expected = np.linalg.eigvalsh(h)
    np.testing.assert_allclose(actual, expected, atol=2e-10, rtol=0,
                               err_msg=label)
    if cg:
        (path / "calc.def").write_text(calc + "CalcType 3\nOutputEigenVec 1\n")
        run(path, "cg", sector, raw_mpi=raw_mpi)
        psi = vector(path, "zvo_eigenvec_0", sector)
        np.testing.assert_allclose(np.linalg.norm(psi), 1, atol=1e-12, rtol=0)
        energy = np.vdot(psi, h @ psi).real
        residual = np.linalg.norm(h @ psi - energy*psi) / max(1, np.linalg.norm(h, 2))
        assert residual <= 1e-10, (label, residual)
        assert abs(energy - expected[0]) <= 1e-10 * max(1, abs(expected[0])), label
    print(f"{label}: {len(words)} eigenvalues" + (" and CG residual passed" if cg else " passed"),
          flush=True)
    return path


LITERAL = pair((0, 0, 1, 1, 1, 1, 2, 0), 1+0j)


def static_checks():
    # The old reader returns [0,0,0] for the first case instead of [-1,0,1].
    for sector in (False, True):
        spectra(f"minimal_{sector}", "Hubbard", [], LITERAL,
                sector=sector, polarized=True, cg=True)

    # Complex contractions on both spins, duplicate rows, reversed partners,
    # and contraction into an already present transfer.
    inter = pair((0, 0, 1, 1, 1, 1, 2, 0), .42+.19j)
    inter += pair((0, 0, 1, 1, 1, 1, 2, 0), .28+.07j)[::-1]
    inter += pair((0, 1, 1, 0, 1, 0, 2, 1), -.31+.23j)[::-1]
    # No middle-index contraction; same-spin products; existing onsite
    # contraction, with explicit chemical potentials also present.
    inter += pair((0, 0, 1, 1, 2, 1, 0, 0), .11+.03j)
    inter += pair((0, 0, 1, 0, 1, 0, 2, 0), -.13+.09j)
    inter += [(0, 0, 1, 1, 1, 1, 0, 0, .37, 0)]
    transfer = [(i, s, i, s, .13*(i+1)+.07*s, 0)
                for i in range(3) for s in range(2)]
    transfer += pair((0, 0, 2, 0), .21-.17j)
    transfer += pair((0, 0, 1, 0), .15+.08j)
    transfer += pair((1, 0, 2, 0), -.12+.05j)
    for model in ("Hubbard", "HubbardGC", "HubbardNConserved",
                  "tJ", "tJGC", "tJNConserved",
                  "Kondo", "KondoGC", "KondoNConserved"):
        # All nine model branches share the reader; fixed-N and projected tJ
        # Hamiltonians are evaluated in the sector implementation.
        if model == "tJGC":
            # TransSym does not support tJGC; its raw projected basis does.
            spectra(model+"_raw", model, transfer, inter, sector=False)
            continue
        spectra(model+"_sector", model, transfer, inter, sector=True, cg=True)
        if model in ("Hubbard", "HubbardGC", "Kondo", "KondoGC"):
            spectra(model+"_raw", model, transfer, inter, sector=False)

    cancel = pair((0, 0, 2, 0), 1+0j)
    for sector in (False, True):
        # Aggregated transfers must retain exact cancellation, including a
        # duplicated explicit input. In the polarized basis H is identically 0.
        spectra(f"cancel_{sector}", "Hubbard",
                [(*row[:-2], row[-2]/2, row[-1]/2) for row in cancel]*2,
                LITERAL, sector=sector, polarized=True)


def schedule(path, times, lists):
    path.write_text("====\nNTimeSteps {}\n====\n====\n====\n".format(len(times))
                   + "".join(f"{t} {len(rows)}\n"
                             + "".join(" ".join(map(str, row))+"\n" for row in rows)
                             for t, rows in zip(times, lists)))


def raw_mpi_check():
    ranks = next((int(v) for k, v in zip(MPI, MPI[1:])
                  if k in ("-np", "-n") and v.isdigit()), 1)
    if ranks not in (4, 16):
        return
    transfer = [(i, s, i, s, .11*(i+1)+.07*s, 0)
                for i in range(4) for s in range(2)]
    transfer += pair((0, 0, 3, 0), .3+.1j)
    transfer += pair((1, 1, 3, 1), -.2+.08j)
    inter = pair((0, 0, 1, 1, 1, 1, 3, 0), .7+.26j)
    spectra("raw_site_mpi", "HubbardGC", transfer, inter,
            sector=False, cg=True, nsite=4, raw_mpi=True)


def time_checks():
    times = [0., .03, .09, .12, .17]
    amplitudes = [0., .21, -.17, 0., .31]
    static_trans = [(i, 0, i, 0, .13*(i+1), 0) for i in range(3)]
    static_trans += pair((0, 0, 1, 0), .3+.1j)
    static_trans += pair((1, 0, 2, 0), .2-.15j)
    # Include static and dynamic contractions together to check the maximum
    # EDTransfer capacity, and switch dynamic terms off after an active slice.
    static_inter = pair((0, 0, 1, 1, 1, 1, 2, 0), .11+.02j)
    dynamic_inter = pair((0, 0, 1, 1, 1, 1, 2, 0), .7+.26j)
    dynamic_inter += pair((0, 0, 0, 1, 0, 1, 2, 0), -.2+.08j)
    dynamic_inter += [(0, 0, 1, 1, 1, 1, 0, 0, .37, 0)]
    dynamic_trans = pair((0, 0, 2, 0), .2-.1j)*2
    dynamic_trans += [(1, 0, 1, 0, .17, 0)]
    for sector in (False, True):
        for family in ("TETwoBody", "TEOneBody"):
            label = f"te_sector{sector}_{family}"
            path, words, calc, mod, names = prepare(
                label, "Hubbard", static_trans, static_inter,
                sector=sector, polarized=True)
            (path / "calc.def").write_text(calc + "CalcType 3\nOutputEigenVec 1\n")
            run(path, "seed", sector)
            psi = vector(path, "zvo_eigenvec_0", sector)
            for file in (path / "output").glob("zvo_eigenvec_0_rank_*.dat"):
                shutil.copyfile(file, file.with_name(file.name.replace("zvo_", "seed_")))
            (path / "calc.def").write_text(
                calc + "CalcType 4\nInputEigenVec 1\nOutputEigenVec 1\n")
            (path / "mod.def").write_text(
                mod + "Lanczos_max 5\nExpandCoef 8\nOutputInterval 1\n"
                "ExpecInterval 1\nTimeSlice .035\nTinit 0\n")
            rows = dynamic_inter if family == "TETwoBody" else dynamic_trans
            lists = [[(*r[:-2], amp*r[-2], amp*r[-1]) for r in rows]
                     if amp else [] for amp in amplitudes]
            # TETwoBody supplies all dynamic one-body capacity in that case.
            # Separate files remain mutually exclusive, as before this fix.
            schedule(path / (family+".def"), times, lists)
            names += f"{family} {family}.def\n"
            (path / "namelist.def").write_text(names + "SpectrumVec seed_eigenvec_0\n")
            continuation = None
            if family == "TETwoBody":
                continuation = path.with_name(label + "_continued")
                shutil.copytree(path, continuation)
            run(path, "evolve", sector)
            static_h = matrix(words, static_trans, static_inter)
            previous_h = static_h
            energies = []
            for step, time in enumerate(times):
                h = static_h + matrix(words,
                                      lists[step] if family == "TEOneBody" else [],
                                      lists[step] if family == "TETwoBody" else [])
                dt = time-times[step-1] if step else 0
                # The sector driver uses the current H for every Taylor
                # power. The existing raw driver caches H(previous)*psi for
                # the first power, then uses the current H. Keep that
                # integration convention separate from the parser oracle.
                first_h = h if sector else previous_h
                term = -1j*dt * (first_h @ psi)
                evolved = psi + term
                for order in range(2, 9):
                    term = (-1j*dt/order) * (h @ term)
                    evolved += term
                psi = evolved / np.linalg.norm(evolved)
                previous_h = h
                energies.append(np.vdot(psi, h @ psi).real)
                # The raw driver overwrites label 0 with its final vector.
                # Step 0 is still checked through its reported energy below.
                if not sector and step == 0:
                    continue
                stem = f"zvo_eigenvec_{step}"
                actual = vector(path, stem, sector)
                overlap = np.vdot(psi, actual)
                np.testing.assert_allclose(actual, psi*overlap/abs(overlap),
                                           atol=2e-10, rtol=0, err_msg=f"{label} step {step}")
            series = np.loadtxt(path / "output/zvo_SS.dat")
            np.testing.assert_allclose(series[:, 0], times, atol=1e-14, rtol=0)
            np.testing.assert_allclose(series[:, 1], energies, atol=2e-10, rtol=0)
            print(f"{label}: step vectors and energies match literal-string evolution", flush=True)
            if continuation is not None:
                # Stop after an active slice; resuming must restore BOTH its
                # four-body term and the one-body contraction before H*psi.
                (continuation / "calc.def").write_text(
                    calc + "CalcType 4\nInputEigenVec 1\nOutputEigenVec 1\nReStart 1\n")
                (continuation / "mod.def").write_text(
                    (path / "mod.def").read_text().replace("Lanczos_max 5", "Lanczos_max 3"))
                run(continuation, "prefix", sector)
                (continuation / "calc.def").write_text(
                    calc + "CalcType 4\nInputEigenVec 0\nOutputEigenVec 1\nReStart 3\n")
                (continuation / "mod.def").write_text((path / "mod.def").read_text())
                run(continuation, "resume", sector)
                for step in (3, 4):
                    np.testing.assert_allclose(
                        vector(continuation, f"zvo_eigenvec_{step}", sector),
                        vector(path, f"zvo_eigenvec_{step}", sector), atol=2e-12, rtol=0)
                print(f"{label}: restart after active slice passed", flush=True)


if __name__ == "__main__":
    print(f"Artifacts: {ROOT}", flush=True)
    static_checks()
    raw_mpi_check()
    time_checks()

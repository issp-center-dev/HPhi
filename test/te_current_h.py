"""Independent four-state oracle for raw TE, including the first Laser step.

Run directly in serial: this checks the integrator, not site decomposition.
te_contract, restart_te_contract and symmetry_te_dynamic cover MPI separately.
The initial state is an arbitrary complex superposition, not an eigenstate.
"""
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

import numpy as np


HPHI = str(Path(sys.argv[1]).resolve())
ROOT = Path(tempfile.mkdtemp(prefix="te_current_h_", dir=".")).resolve()
ENV = dict(os.environ, OMP_NUM_THREADS="1")
# Two sites, one up and one down electron, in ascending Fock-word order.
WORDS = [3, 6, 9, 12]
INITIAL = np.array([1, .2+.3j, -.4+.2j, .6j])
INITIAL /= np.linalg.norm(INITIAL)
TOL = 2e-10


def pair(indices, value):
    return [(*indices, value.real, value.imag),
            (*indices[2:], *indices[:2], value.real, -value.imag)]


TRANS = [(i, s, i, s, .17*(i+1)+.09*s, 0)
         for i in range(2) for s in range(2)]
TRANS += sum((pair((0, s, 1, s), .7+.15j) for s in range(2)), [])
INTER = [(i, 0, i, 0, i, 1, i, 1, .8+.5*i, 0) for i in range(2)]
DYNAMIC_TRANS = pair((0, 0, 1, 0), -.3+.4j) + [(0, 1, 0, 1, .6, 0)]
DYNAMIC_INTER = [(0, 0, 0, 0, 0, 1, 0, 1, 1.7, 0)]


def matrix(transfer, interactions):
    """Apply literal fermionic operator strings, without HPhi normalization."""
    h = np.zeros((4, 4), complex)
    for rows, prefactor in ((transfer, -1), (interactions, 1)):
        for row in rows:
            indices = row[:-2]
            operators = [(j % 2 == 0, 2*indices[2*j]+indices[2*j+1])
                         for j in range(len(indices)//2)]
            for col, word in enumerate(WORDS):
                value = prefactor * complex(*row[-2:])
                for create, orbital in reversed(operators):
                    if bool(word & (1 << orbital)) == create:
                        break
                    value *= (-1)**((word & ((1 << orbital)-1)).bit_count())
                    word ^= 1 << orbital
                else:
                    if word in WORDS:
                        h[WORDS.index(word), col] += value
    np.testing.assert_allclose(h, h.conj().T, atol=1e-14, rtol=0)
    return h


def scaled(rows, scale):
    return [(*r[:-2], scale*r[-2], scale*r[-1]) for r in rows] if scale else []


def definition(path, rows, label="NItems", count=None):
    count = len(rows) if count is None else count
    path.write_text(f"====\n{label} {count}\n====\n====\n====\n" +
                    "".join(" ".join(map(str, r))+"\n" for r in rows))


def run(path, label):
    log = path / (label+".log")
    with log.open("w") as fp:
        result = subprocess.run([HPHI, "-e", "namelist.def"], cwd=path,
                                env=ENV, stdout=fp, stderr=subprocess.STDOUT,
                                timeout=60)
    assert result.returncode == 0, (path, log.read_text())


def vector(path, step=0):
    data = (path / f"output/zvo_eigenvec_{step}_rank_0.dat").read_bytes()
    _, size = struct.unpack_from("=iQ", data)
    assert size == 4 and len(data) == 12+16*(size+1)
    return np.frombuffer(data, dtype="=c16", offset=12)[1:].copy()


def state_error(expected, actual):
    overlap = np.vdot(expected, actual)
    assert abs(overlap) > 0
    return np.linalg.norm(actual - expected * overlap/abs(overlap))


def exact_step(h, dt, psi, order=None):
    # Diagonalization is independent of the production matvec recurrence.
    values, vectors = np.linalg.eigh(h)
    z = -1j*dt*values
    factor = (np.exp(z) if order is None else
              sum(z**k / math.factorial(k) for k in range(order+1)))
    evolved = vectors @ (factor * (vectors.conj().T @ psi))
    norm = np.linalg.norm(evolved)
    return evolved/norm, norm


def prepare(label, family, times, amplitudes, order=12, smooth=False):
    path = ROOT / label
    path.mkdir()
    (path / "output").mkdir()
    (path / "output/seed_rank_0.dat").write_bytes(
        struct.pack("=iQ", 0, 4) + np.r_[0j, INITIAL].astype("=c16").tobytes())
    (path / "calc.def").write_text(
        "CalcModel 0\nCalcType 4\nOutputMode 0\nOutputDataHead 1\n"
        "InputEigenVec 1\nOutputEigenVec 1\nReStart 0\n")
    (path / "mod.def").write_text(
        "====\nModel_Parameters 0\n====\n====\n====\n"
        "CDataFileHead zvo\nCParaFileHead zqp\n====\n"
        "Nsite 2\nNcond 2\n2Sz 0\ninitial_iv -1\nexct 1\n"
        "LanczosEps 14\nLanczosTarget 1\nLargeValue 100\n"
        f"Lanczos_max {len(times)}\nExpandCoef {order}\n"
        "OutputInterval 1\nExpecInterval 1\nTimeSlice .06\nTinit .13\n")
    definition(path / "loc.def", [(0, 0), (1, 0)], "NlocalSpin", count=0)
    inter = [] if smooth else INTER
    definition(path / "trans.def", TRANS)
    definition(path / "inter.def", inter)
    names = ("CalcMod calc.def\nModPara mod.def\nLocSpin loc.def\n"
             "Trans trans.def\nInterAll inter.def\nSpectrumVec seed\n")
    drive = TRANS if smooth else DYNAMIC_TRANS if family == "one" else DYNAMIC_INTER
    lists = [scaled(drive, amp) if family in ("one", "two") else []
             for amp in amplitudes]
    table = f"====\nNTimeSteps {len(times)}\n====\n====\n====\n"
    for time, rows in zip(times, lists):
        table += f"{time:.17g} {len(rows)}\n"
        table += "".join(" ".join(map(str, r))+"\n" for r in rows)
    (path / "times.def").write_text(table)
    names += ("TETwoBody" if family == "two" else "TEOneBody") + " times.def\n"
    if family == "laser":
        definition(path / "laser.def", list(enumerate([2, .8, 1, 1, 0, 2, 1, 1, 0])), "NLaser")
        names += "Laser laser.def\n"
    (path / "namelist.def").write_text(names)
    matrices = []
    for step, rows in enumerate(lists):
        trans = TRANS
        if family == "laser":
            # Raw Laser samples H at the left end and outputs the right end.
            a = .8 * (.13 + step*.06)
            trans = [(*r[:-2], (complex(*r[-2:])*np.exp(-1j*a*(r[0]-r[2]))).real,
                      (complex(*r[-2:])*np.exp(-1j*a*(r[0]-r[2]))).imag) for r in TRANS]
        matrices.append(matrix(trans + (rows if family == "one" else []),
                               inter + (rows if family == "two" else [])))
    return path, matrices


def check(label, family, times, amplitudes, order=12, smooth=False):
    path, matrices = prepare(label, family, times, amplitudes, order, smooth)
    run(path, "evolve")
    psi = INITIAL.copy()
    energies, norms, states = [], [], []
    for step, h in enumerate(matrices):
        dt = .06 if family == "laser" else times[step]-times[step-1] if step else 0
        # Order 12: compare with exp, not the same Taylor implementation.
        # Order 1: independently evaluate the polynomial in the eigenbasis.
        psi, norm = exact_step(h, dt, psi, None if order == 12 else order)
        states.append(psi)
        energies.append(np.vdot(psi, h @ psi).real)
        norms.append(norm)
    errors = [state_error(states[i], vector(path, i)) for i in range(1, len(times))]
    errors.append(state_error(psi, vector(path)))  # raw label 0 is the final alias
    energy = np.loadtxt(path / "output/zvo_SS.dat", ndmin=2)
    norm = np.loadtxt(path / "output/zvo_Norm.dat", ndmin=2)
    labels = .13 + (np.arange(len(times))+1)*.06 if family == "laser" else times
    np.testing.assert_allclose(energy[:, 0], labels, atol=1e-14, rtol=0)
    energy_error = np.max(abs(energy[:, 1]-energies))
    print(f"{label}: state={max(errors):.6e}, energy={energy_error:.6e}, "
          f"norm={np.max(abs(norm[:, 1]-norms)):.6e}", flush=True)
    assert max(errors) < TOL, (label, errors)
    np.testing.assert_allclose(energy[:, 1], energies, atol=TOL, rtol=0)
    np.testing.assert_allclose(norm[:, 1], norms, atol=TOL, rtol=0)
    return path


def restart_checks(reference):
    for v3 in (False, True):
        for split in (1, 3):
            path = reference.with_name(reference.name+f"_v3{v3}_split{split}")
            shutil.copytree(reference, path)
            calc = (path / "calc.def").read_text()
            mod = (path / "mod.def").read_text()
            assert "Lanczos_max 5\n" in mod and "ReStart 0\n" in calc
            (path / "calc.def").write_text(calc.replace("ReStart 0", f"ReStart {int(v3)}"))
            (path / "mod.def").write_text(mod.replace("Lanczos_max 5", f"Lanczos_max {split}"))
            run(path, "prefix")
            if not v3:
                shutil.copyfile(path / "output/zvo_eigenvec_0_rank_0.dat", path / "output/seed_rank_0.dat")
            (path / "calc.def").write_text(calc.replace("ReStart 0", "ReStart 3"))
            (path / "mod.def").write_text(mod)
            run(path, "continued")
            np.testing.assert_allclose(vector(path), vector(reference), atol=2e-12, rtol=0)
            for name in ("SS", "Norm", "Flct"):
                np.testing.assert_allclose(np.loadtxt(path / f"output/zvo_{name}.dat"),
                                           np.loadtxt(reference / f"output/zvo_{name}.dat"),
                                           atol=2e-12, rtol=0)
    print(f"{reference.name}: first/middle legacy and v3 restart passed", flush=True)


def convergence():
    errors = []
    # H(t)=(1+t)K commutes at all times: integral_0^T (1+t)dt=T+T^2/2.
    total = .4
    exact, _ = exact_step(matrix(TRANS, []), total+total**2/2, INITIAL)
    for n in (4, 8, 16):
        times = np.linspace(0, total, n+1)
        path = check(f"smooth_n{n}", "one", times, times, smooth=True)
        error = state_error(exact, vector(path))
        errors.append(error)
        print(f"smooth: dt={total/n:.6g}, continuous_state_error={error:.9e}", flush=True)
    ratios = np.array(errors[:-1])/errors[1:]
    assert np.all((ratios > 1.95) & (ratios < 2.05)), ratios


def main():
    times = [0, .07, .12, .12, .24]
    amps = [0, .9, -.5, 0, .8]
    cases = [("static", "static", times, amps, 12),
             ("zero_dt", "one", [0], [.8], 12),
             ("one", "one", times, amps, 12),
             ("two", "two", times, amps, 12),
             ("laser_first", "laser", [0], [0], 12),
             ("laser", "laser", times, amps, 12),
             ("one_order1", "one", times, amps, 1),
             ("laser_order1", "laser", times, amps, 1)]
    failures = []
    for label, *args in cases:
        try:
            path = check(label, *args)
            if label in ("static", "one", "two", "laser"):
                restart_checks(path)
        except (AssertionError, subprocess.TimeoutExpired) as error:
            failures.append((label, str(error)))
    try:
        convergence()
    except (AssertionError, subprocess.TimeoutExpired) as error:
        failures.append(("convergence", str(error)))
    assert not failures, (ROOT, failures)
    print("raw TE current-H oracle passed", flush=True)


if __name__ == "__main__":
    main()

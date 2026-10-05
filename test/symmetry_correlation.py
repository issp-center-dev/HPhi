"""Sector correlation functions against an independent physical-basis reference.

The reference applies operator strings ket by ket in the fixed-quantum-number
physical basis (final-state projection for tJ). Its sign convention is
cross-checked against the Jordan-Wigner tensor products of
symmetry_general_terms for every model at a size where the full space is small.
"""
import cmath
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import symmetry_general_terms as fixture

fixture.ROOT = Path("symmetry_correlation")
fixture.ROOT.mkdir(exist_ok=True)
ATOL = 1e-8

CASES = {
    "Spin": dict(calc=1, L=6, nup=3, ndown=3, quantum="2Sz 0\n", group="dihedral"),
    "SpinlessFermion": dict(calc=7, L=6, nup=3, ndown=0, quantum="Ncond 3\n", group="translation"),
    "Hubbard": dict(calc=0, L=4, nup=2, ndown=2, quantum="Nup 2\nNdown 2\n", group="translation"),
    "tJ": dict(calc=9, L=6, nup=2, ndown=2, quantum="Nup 2\nNdown 2\n", group="translation"),
}
FILES = {"one": "cisajs", "two": "cisajscktalt", "three": "ThreeBody", "four": "FourBody",
         "six": "SixBody", "nbody": "NBodyG"}
KEYWORDS = {"one": "OneBodyG", "two": "TwoBodyG", "three": "ThreeBodyG", "four": "FourBodyG",
            "six": "SixBodyG", "nbody": "NBodyG"}


def width_of(model):
    return 2 if model in ("Hubbard", "tJ") else 1


def physical_states(model, L, nup, ndown):
    if model in ("Spin", "SpinlessFermion"):
        return [s for s in range(1 << L) if bin(s).count("1") == nup]
    states = [s for s in range(1 << (2 * L))
              if sum((s >> (2 * i)) & 1 for i in range(L)) == nup
              and sum((s >> (2 * i + 1)) & 1 for i in range(L)) == ndown]
    if model == "tJ":
        states = [s for s in states if all(((s >> (2 * i)) & 3) != 3 for i in range(L))]
    return states


def apply_factors(model, L, state, factors):
    """Apply c^dag_{i s} c_{j t} factors right to left; (state, sign) or None."""
    width = width_of(model)
    sign = 1
    for i, s, j, t in reversed(factors):
        if model == "Spin":
            if i != j or (state >> i) & 1 != t:
                return None
            state = (state & ~(1 << i)) | (s << i)
            continue
        for orbital, create in ((width * j + t, False), (width * i + s, True)):
            mask = 1 << orbital
            if bool(state & mask) == create:
                return None
            sign *= (-1) ** bin(state & (mask - 1)).count("1")
            state ^= mask
    if model == "tJ" and any(((state >> (2 * k)) & 3) == 3 for k in range(L)):
        return None
    return state, sign


def matrix(model, L, states, factors, coefficient=1.0):
    index = {s: n for n, s in enumerate(states)}
    out = np.zeros((len(states), len(states)), dtype=complex)
    for col, s in enumerate(states):
        result = apply_factors(model, L, s, factors)
        if result is None:
            continue
        row = index.get(result[0])
        if row is not None:
            out[row, col] += coefficient * result[1]
    return out


def crosscheck_builder(model, L, nup, ndown):
    """apply_factors must agree with the Jordan-Wigner tensor products of the fixture."""
    width = width_of(model)
    full = width * L
    states = physical_states(model, L, nup, ndown)
    if model == "Spin":
        def factor(i, s, j, t):
            unit = np.zeros((2, 2))
            unit[s, t] = 1
            return fixture.tensor_operator(full, i, unit)
    else:
        annihilator = {o: fixture.tensor_operator(full, o, np.array([[0, 1], [0, 0]]), True) for o in range(full)}

        def factor(i, s, j, t):
            return annihilator[width * i + s].T @ annihilator[width * j + t]
    rng = np.random.default_rng(7)
    for _ in range(40):
        factors = []
        for _ in range(int(rng.integers(1, 4))):
            i = int(rng.integers(L))
            j = i if model == "Spin" else int(rng.integers(L))
            s, t = (0, 0) if model == "SpinlessFermion" else (int(rng.integers(2)), int(rng.integers(2)))
            factors.append((i, s, j, t))
        dense = np.eye(1 << full)
        for f in factors:
            dense = dense @ factor(*f)
        np.testing.assert_allclose(matrix(model, L, states, factors), dense[np.ix_(states, states)], atol=1e-12)


def hamiltonian(model, L):
    """Expert-mode families; dense meaning per family as pinned by symmetry_general_terms."""
    fam = {}

    def add(name, row):
        fam.setdefault(name, []).append(row)
    for i in range(L):
        j = (i + 1) % L
        if model == "Spin":
            add("Exchange", [i, j, 0.5])
            for s in (0, 1):
                for t in (0, 1):
                    add("InterAll", [i, s, i, s, j, t, j, t, 0.25 if s == t else -0.25, 0])
        elif model == "SpinlessFermion":
            add("Trans", [i, 0, j, 0, 1.0, 0])
            add("Trans", [j, 0, i, 0, 1.0, 0])
            add("InterAll", [i, 0, i, 0, j, 0, j, 0, 1.0, 0])
        else:
            for s in (0, 1):
                add("Trans", [i, s, j, s, 1.0, 0])
                add("Trans", [j, s, i, s, 1.0, 0])
            if model == "Hubbard":
                add("CoulombIntra", [i, 4.0])
            else:
                add("InterAll", [i, 0, i, 1, j, 1, j, 0, 0.25, 0])
                add("InterAll", [j, 0, j, 1, i, 1, i, 0, 0.25, 0])
                add("InterAll", [i, 0, i, 0, j, 1, j, 1, -0.25, 0])
                add("InterAll", [i, 1, i, 1, j, 0, j, 0, -0.25, 0])
    return fam


def dense_hamiltonian(model, L, states, fam):
    H = np.zeros((len(states),) * 2, dtype=complex)
    for name, rows in fam.items():
        for row in rows:
            if name == "Trans":
                i, s, j, t, re, im = row
                H -= complex(re, im) * matrix(model, L, states, [(i, s, j, t)])
            elif name == "InterAll":
                i, s, j, t, k, u, l, v, re, im = row
                H += complex(re, im) * matrix(model, L, states, [(i, s, j, t), (k, u, l, v)])
            elif name == "Exchange":
                i, j, J = row
                H += J * (matrix(model, L, states, [(i, 0, i, 1), (j, 1, j, 0)])
                          + matrix(model, L, states, [(j, 0, j, 1), (i, 1, i, 0)]))
            elif name == "CoulombIntra":
                i, U = row
                H += U * matrix(model, L, states, [(i, 0, i, 0), (i, 1, i, 1)])
    np.testing.assert_allclose(H, H.conj().T, atol=1e-12)
    return H


def operations(L, group):
    ops = []
    for p in ((0,) if group == "translation" else (0, 1)):
        for g in range(L):
            perm = [((L - 1 - i) if p else i) for i in range(L)]
            ops.append(([(x + g) % L for x in perm], g, p))
    return ops


def transform(model, L, state, perm):
    width = width_of(model)
    occupied = [o for o in range(width * L) if (state >> o) & 1]
    moved = [width * perm[o // width] + o % width for o in occupied]
    sign = 1 if model == "Spin" else (-1) ** sum(a > b for n, a in enumerate(moved) for b in moved[n + 1:])
    return sum(1 << o for o in moved), sign


def characters(L, ops, k, parity):
    return [cmath.exp(-2j * cmath.pi * k * g / L) * (parity ** p) for _, g, p in ops]


def projector(model, L, states, ops, chars):
    index = {s: n for n, s in enumerate(states)}
    P = np.zeros((len(states),) * 2, dtype=complex)
    for (perm, _, _), ch in zip(ops, chars):
        for col, s in enumerate(states):
            t, sign = transform(model, L, s, perm)
            P[index[t], col] += ch.conjugate() * sign / len(ops)
    return P


def sector_basis(P):
    weights, vectors = np.linalg.eigh(P)
    return vectors[:, weights > 0.5]


def sectors(L, group):
    if group == "translation":
        return [(k, 1) for k in range(L)]
    return [(0, 1), (0, -1), (L // 2, 1), (L // 2, -1)]


def factors_of(row, factors):
    return [tuple(row[4 * f:4 * f + 4]) for f in range(factors)]


def green_requests(model, L):
    spins = (0,) if model == "SpinlessFermion" else (0, 1)
    one = [[i, s, j, t] for i in range(L) for j in range(L) for s in spins for t in spins
           if model != "Spin" or i == j]
    two = [[i, s, i, s, j, t, j, t] for i in range(L) for j in range(L) for s in spins for t in spins]
    if len(spins) == 2:
        two += [[i, 0, i, 1, j, 1, j, 0] for i in range(L) for j in range(L) if i != j]
    if model == "Hubbard":
        two += [[0, 0, 1, 0, 0, 1, 1, 1], [1, 0, 2, 0, 1, 1, 2, 1]]
    if model == "SpinlessFermion":
        two += [[0, 0, 1, 0, 2, 0, 2, 0], [0, 0, 2, 0, 1, 0, 1, 0]]
    if model == "Spin":
        nbody = [[0, 0, 0, 1, 1, 1, 1, 0, 2, 0, 2, 0], [0, 1, 0, 1, 1, 0, 1, 0, 2, 1, 2, 1],
                 [1, 0, 1, 1, 2, 1, 2, 0, 0, 0, 0, 0]]
    elif model == "SpinlessFermion":
        nbody = [[0, 0, 1, 0, 1, 0, 2, 0, 2, 0, 0, 0], [0, 0, 0, 0, 1, 0, 1, 0, 2, 0, 2, 0],
                 [0, 0, 2, 0, 3, 0, 3, 0, 1, 0, 0, 0]]
    else:
        nbody = [[0, 0, 1, 0, 1, 1, 1, 1, 2, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1, 2, 0, 2, 0],
                 [0, 0, 0, 1, 1, 1, 1, 0, 2, 0, 2, 0]]
    requests = {"one": one, "two": two, "nbody": [[3] + row for row in nbody]}
    if model == "Hubbard":
        requests["three"] = [[0, 0, 1, 0, 1, 1, 1, 1, 2, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1, 2, 0, 2, 0]]
        requests["four"] = [[0, 0, 1, 0, 1, 0, 0, 0, 2, 1, 3, 1, 3, 1, 2, 1]]
        # hopping pairs and repeated densities: nonzero for Nup = Ndown = 2 (asserted below)
        requests["six"] = [[0, 0, 1, 0, 1, 0, 0, 0, 2, 1, 3, 1, 3, 1, 2, 1, 0, 0, 0, 0, 1, 1, 1, 1]]
    return requests


def request_factors(kind, row):
    if kind == "nbody":
        return factors_of(row[1:], row[0])
    return factors_of(row, {"one": 1, "two": 2, "three": 3, "four": 4, "six": 6}[kind])


def expected_values(model, L, states, gs, kind, rows):
    return np.array([np.vdot(gs, matrix(model, L, states, request_factors(kind, row)) @ gs) for row in rows])


def write_inputs(path, model, case, fam, ops, chars, method, requests, aggregate=0):
    L = case["L"]
    fixture.definition(path, "loc.def", [[i, int(model == "Spin")] for i in range(L)], L if model == "Spin" else 0)
    namelist = "CalcMod calc.def\nModPara mod.def\nLocSpin loc.def\n"
    for family, rows in fam.items():
        fixture.definition(path, family + ".def", rows)
        namelist += "{} {}.def\n".format(family, family)
    group = [[o, ch.real, ch.imag] for o, ch in enumerate(chars)]
    for o, (perm, _, _) in enumerate(ops):
        group.extend([[o, i, perm[i], 1] for i in range(L)])
    fixture.definition(path, "group.def", group, len(ops))
    (path / "group.def").write_text((path / "group.def").read_text().replace("NItems", "NQPTrans"))
    for kind, rows in requests.items():
        fixture.definition(path, FILES[kind] + ".def", rows)
        namelist += "{} {}.def\n".format(KEYWORDS[kind], FILES[kind])
    (path / "raw.def").write_text(namelist)
    (path / "sym.def").write_text(namelist + "TransSym group.def\n")
    (path / "calc.def").write_text("CalcType {}\nCalcModel {}\nOutputMode 0\nOutputDataHead 1\nOutputGreenFormat {}\n".format(
        method, case["calc"], aggregate))
    (path / "mod.def").write_text(
        "====\nModel_Parameters 0\n====\n====\n====\nCDataFileHead zvo\nCParaFileHead zqp\n====\n"
        "Nsite {}\n{}Lanczos_max 2000\ninitial_iv -1\nexct 1\nLanczosEps 20\nLanczosTarget 1\nLargeValue 100\nPreCG 0\n".format(
            L, case["quantum"]))


def output_name(method, kind, aggregate=False):
    if aggregate:
        candidates = sorted(Path("output").glob("zvo_{}_*.dat".format(FILES[kind])))
        return candidates
    return "zvo_{}_eigen0.dat".format(FILES[kind]) if method == 3 else "zvo_{}.dat".format(FILES[kind])


def check_file(path, name, rows, expected, prefix=0, label=""):
    data = np.loadtxt(path / "output" / name, ndmin=2)
    assert data.shape == (len(rows), prefix + len(rows[0]) + 2), (label, name, data.shape)
    np.testing.assert_array_equal(data[:, prefix:-2], np.array(rows, dtype=float), err_msg=label + name)
    values = data[:, -2] + 1j * data[:, -1]
    np.testing.assert_allclose(values, expected, atol=ATOL, rtol=0, err_msg=label + name)


def check_fixed_body_spacing(path, name, model, kind, rows):
    factors = {"two": 2, "three": 3, "four": 4, "six": 6}[kind]
    lines = (path / "output" / name).read_text().splitlines()
    assert len(lines) == len(rows), (model, kind, len(lines), len(rows))
    for line, row in zip(lines, rows):
        value_fields = line.split()[-2:]
        index_fields = ""
        for k, value in enumerate(row):
            if (factors >= 4 and k == 12) or (factors >= 6 and k == 16):
                index_fields += " "
            index_fields += " {:4d}".format(value)
        trailing = " " if (factors >= 3 or model == "Spin") else ""
        expected = index_fields + " {} {}{}".format(value_fields[0], value_fields[1], trailing)
        assert line == expected, (model, kind, repr(line), repr(expected))


def reset(path):
    if path.exists():
        shutil.rmtree(str(path))
    path.mkdir(parents=True)


def ground_state(model, case, H, states, ops):
    L = case["L"]
    global_min = np.linalg.eigvalsh(H)[0]
    best = None
    for k, parity in sectors(L, case["group"]):
        chars = characters(L, ops, k, parity)
        P = projector(model, L, states, ops, chars)
        np.testing.assert_allclose(H @ P, P @ H, atol=1e-12)
        B = sector_basis(P)
        if B.shape[1] == 0:
            continue
        E, V = np.linalg.eigh(B.conj().T @ H @ B)
        if best is None or E[0] < best["E"][0] - 1e-12:
            best = dict(E=E, k=k, parity=parity, chars=chars, basis=B, vectors=V)
    assert abs(best["E"][0] - global_min) < 1e-10
    assert len(best["E"]) == 1 or best["E"][1] - best["E"][0] > 1e-6, best["E"][:2]
    best["gs"] = best["basis"] @ best["vectors"][:, 0]
    return best


def run_and_check(path, model, case, fam, ops, sector, requests, expected, method, layout, label, aggregate=0):
    write_inputs(path, model, case, fam, ops, sector["chars"], method, requests, aggregate)
    if (path / "output").exists():
        shutil.rmtree(str(path / "output"))
    energies = fixture.run(path, label, layout=layout)
    np.testing.assert_allclose(energies[0], sector["E"][0], atol=3e-8, rtol=0, err_msg=label)
    stats = [line for line in (path / (label + ".log")).read_text().splitlines()
             if line.startswith("Symmetry correlation:")]
    assert stats, label
    multiwave_models = os.environ.get("HPHI_TEST_EXPECT_MULTIWAVE", "").split(",")
    if "*" in multiwave_models or path.name in multiwave_models:
        assert any(int(line.split("waves=")[1].split()[0]) >= 2 for line in stats), (label, stats)
    if os.environ.get("HPHI_TEST_EXPECT_THREADS"):
        wanted = int(os.environ["HPHI_TEST_EXPECT_THREADS"])
        assert all(int(line.split("threads=")[1].split()[0]) >= wanted for line in stats), (label, stats)
    for kind, rows in requests.items():
        if aggregate:
            files = sorted((path / "output").glob("zvo_{}_*.dat".format(FILES[kind])))
            assert len(files) == 1, (label, kind, files)
            check_file(path, files[0].name, rows, expected[kind], prefix=1, label=label + " ")
            data = np.loadtxt(files[0], ndmin=2)
            assert np.all(data[:, 0] == 0), (label, "eigen prefix")
        else:
            name = output_name(method, kind)
            check_file(path, name, rows, expected[kind], label=label + " ")
            if kind in ("two", "three", "four", "six"):
                check_fixed_body_spacing(path, name, model, kind, rows)
    shutil.copytree(str(path / "output"), str(path / label))


def test_model(model):
    case = CASES[model]
    L = case["L"]
    path = fixture.ROOT / model
    reset(path)
    states = physical_states(model, L, case["nup"], case["ndown"])
    crosscheck_builder(model, 4 if model == "tJ" else L, 2 if model == "tJ" else case["nup"],
                       2 if model == "tJ" else case["ndown"])
    fam = hamiltonian(model, L)
    H = dense_hamiltonian(model, L, states, fam)
    ops = operations(L, case["group"])
    sector = ground_state(model, case, H, states, ops)
    if model == "Spin":
        assert sector["k"] == 3 and abs(sector["E"][0] + 2.8027756377319952) < 1e-10, sector["k"]
    requests = green_requests(model, L)
    expected = {kind: expected_values(model, L, states, sector["gs"], kind, rows) for kind, rows in requests.items()}
    # Every family must contain a row whose operator and expectation value are nonzero,
    # so that an implementation writing zeros cannot pass.
    for kind, rows in requests.items():
        norms = [np.linalg.norm(matrix(model, L, states, request_factors(kind, row))) for row in rows]
        assert max(norms) > 1e-12, (model, kind, norms)
        assert np.max(np.abs(expected[kind])) > 1e-6, (model, kind, expected[kind])
    # Rows related by a symmetry operation carry the same value in the reference, too.
    for perm, _, _ in ops:
        for row, value in zip(requests["one"], expected["one"]):
            moved = [perm[row[0]], row[1], perm[row[2]], row[3]]
            if moved in requests["one"]:
                np.testing.assert_allclose(expected["one"][requests["one"].index(moved)], value, atol=1e-12)
    layouts = tuple(filter(None, os.environ.get(
        "HPHI_TEST_LAYOUTS", "replicated,distributed").split(",")))
    assert layouts and all(layout in ("replicated", "distributed") for layout in layouts), layouts
    for layout in layouts:
        run_and_check(path, model, case, fam, ops, sector, requests, expected, 3, layout, "cg_" + layout)
    if os.environ.get("HPHI_TEST_PRIMARY_ONLY"):
        print("{}: primary sector correlation run verified".format(model))
        return
    if model in ("Spin", "SpinlessFermion"):
        run_and_check(path, model, case, fam, ops, sector, requests, expected, 0, "distributed", "lanczos")
    if model == "Hubbard":
        run_and_check(path, model, case, fam, ops, sector, requests, expected, 3, "distributed", "cg_aggregate", aggregate=1)
        higher = {kind: requests[kind] for kind in ("three", "four", "six")}
        run_and_check(path, model, case, fam, ops, sector, higher,
                      {kind: expected[kind] for kind in higher}, 3, "distributed", "cg_higher_only")
    if model == "SpinlessFermion":
        chars = characters(L, ops, 1, 1)
        P = projector(model, L, states, ops, chars)
        B = sector_basis(P)
        assert B.shape[1] == 3, B.shape
        E, V = np.linalg.eigh(B.conj().T @ H @ B)
        assert E[1] - E[0] > 1e-6
        complex_sector = dict(E=E, chars=chars, gs=B @ V[:, 0])
        one = [[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 2, 0]]
        values = expected_values(model, L, states, complex_sector["gs"], "one", one)
        assert np.max(np.abs(values.imag)) > 1e-3, values
        for layout in ("replicated", "distributed"):
            run_and_check(path, model, case, fam, ops, complex_sector, {"one": one}, {"one": values},
                          3, layout, "cg_k1_" + layout)
    if model == "Spin":
        spin_raw_comparison(path, model, case, fam, states, H, requests, expected)
        spin_reader_drop(path, model, case, fam, requests)
    if model == "tJ":
        fermion_raw_comparison(path, model, case, fam, ops, sector, states, requests)
    negative_tests(path, model, case, fam, ops, sector, requests)
    print("{}: sector correlation functions verified".format(model))


def spin_raw_comparison(path, model, case, fam, states, H, requests, expected):
    """Translation-only group at k=pi: the raw ground state lives in this sector."""
    L = case["L"]
    ops = operations(L, "translation")
    chars = characters(L, ops, 3, 1)
    P = projector(model, L, states, ops, chars)
    B = sector_basis(P)
    E, V = np.linalg.eigh(B.conj().T @ H @ B)
    assert E[1] - E[0] > 1e-6
    sector = dict(E=E, chars=chars, gs=B @ V[:, 0])
    subset = {"one": requests["one"], "two": requests["two"]}
    values = {kind: expected_values(model, L, states, sector["gs"], kind, rows) for kind, rows in subset.items()}
    run_and_check(path, model, case, fam, ops, sector, subset, values, 3, "distributed", "cg_translation_kpi")
    write_inputs(path, model, case, fam, ops, chars, 3, subset)
    shutil.rmtree(str(path / "output"))
    raw = fixture.run(path, "raw_kpi", symmetry=False, layout=None)
    np.testing.assert_allclose(raw[0], E[0], atol=3e-8, rtol=0)
    for kind in subset:
        name = output_name(3, kind)
        raw_text = (path / "output" / name).read_text().split()
        sector_text = (path / "cg_translation_kpi" / name).read_text().split()
        assert len(raw_text) == len(sector_text), name
        np.testing.assert_allclose(np.array(raw_text, dtype=float), np.array(sector_text, dtype=float),
                                   atol=1e-8, rtol=0, err_msg="raw vs sector " + name)


def spin_reader_drop(path, model, case, fam, requests):
    """Non-on-site Spin OneBodyG rows are dropped by the reader in both raw and sector runs."""
    L = case["L"]
    ops = operations(L, "translation")
    chars = characters(L, ops, 3, 1)
    rows = requests["one"] + [[0, 0, 1, 0], [2, 1, 3, 1]]
    write_inputs(path, model, case, fam, ops, chars, 3, {"one": rows})
    outputs = {}
    for label, symmetry in (("drop_raw", False), ("drop_sector", True)):
        if (path / "output").exists():
            shutil.rmtree(str(path / "output"))
        fixture.run(path, label, symmetry=symmetry, layout=None if not symmetry else "distributed")
        data = np.loadtxt(path / "output" / output_name(3, "one"), ndmin=2)
        assert data.shape[0] == len(requests["one"]), (label, data.shape)
        np.testing.assert_array_equal(data[:, :4], np.array(requests["one"], dtype=float))
        outputs[label] = data
    np.testing.assert_allclose(outputs["drop_raw"], outputs["drop_sector"], atol=1e-8, rtol=0)


def fermion_raw_comparison(path, model, case, fam, ops, sector, states, requests):
    """Check the projected tJ convention against a raw-basis calculation."""
    def conserves_spin(kind, row):
        factors = request_factors(kind, row)
        return sorted(factor[1] for factor in factors) == sorted(factor[3] for factor in factors)

    subset = {kind: [row for row in requests[kind] if conserves_spin(kind, row)]
              for kind in ("one", "two", "nbody")}
    assert all(subset.values())
    expected = {kind: expected_values(model, case["L"], states, sector["gs"], kind, rows)
                for kind, rows in subset.items()}
    run_and_check(path, model, case, fam, ops, sector, subset, expected, 3,
                  "distributed", "cg_raw_reference")
    write_inputs(path, model, case, fam, ops, sector["chars"], 3, subset)
    if (path / "output").exists():
        shutil.rmtree(str(path / "output"))
    raw = fixture.run(path, "raw_tj", symmetry=False, layout=None)
    np.testing.assert_allclose(raw[0], sector["E"][0], atol=3e-8, rtol=0)
    for kind, rows in subset.items():
        name = output_name(3, kind)
        raw_data = np.loadtxt(path / "output" / name, ndmin=2)
        sector_data = np.loadtxt(path / "cg_raw_reference" / name, ndmin=2)
        assert raw_data.shape == sector_data.shape == (len(rows), len(rows[0]) + 2), (kind, raw_data.shape)
        np.testing.assert_array_equal(raw_data[:, :-2], sector_data[:, :-2])
        np.testing.assert_allclose(raw_data[:, -2:], sector_data[:, -2:], atol=ATOL, rtol=0,
                                   err_msg="raw fermion vs sector " + name)
    shutil.copytree(str(path / "output"), str(path / "raw_tj"))


def negative_tests(path, model, case, fam, ops, sector, requests):
    L = case["L"]
    write_inputs(path, model, case, fam, ops, sector["chars"], 2, {"one": requests["one"]})
    fixture.run(path, "reject_fulldiag", fail="sector FullDiag outputs eigenvalues only")
    if model == "Spin":
        write_inputs(path, model, case, fam, ops, sector["chars"], 3,
                     {"one": requests["one"], "three": [[0, 0, 0, 1, 1, 1, 1, 0, 2, 0, 2, 0]]})
        fixture.run(path, "reject_threebody", fail="ThreeBodyG/FourBodyG/SixBodyG for canonical Spin is not supported")
        # The reader itself rejects non-on-site Spin TwoBodyG rows (CheckFormatForSpinInt, exitMPI).
        write_inputs(path, model, case, fam, ops, sector["chars"], 3, {"two": [[0, 0, 1, 0, 2, 0, 2, 0]]})
        fixture.run(path, "reject_spin_twobody_format", fail="i=j and k=l must be satisfied")
    if model == "Hubbard":
        write_inputs(path, model, case, fam, ops, sector["chars"], 3, {"one": requests["one"]})
        fixture.definition(path, "anomalous.def", [[0, 0, 0, 1, 1]])
        (path / "sym.def").write_text((path / "sym.def").read_text() + "AnomalousG anomalous.def\n")
        fixture.run(path, "reject_anomalous", fail="AnomalousG")


if __name__ == "__main__":
    models = tuple(filter(None, os.environ.get(
        "HPHI_TEST_MODELS", "Spin,SpinlessFermion,Hubbard,tJ").split(",")))
    assert models and all(model in CASES for model in models), models
    for model in models:
        test_model(model)

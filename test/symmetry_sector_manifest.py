"""Independent orbit enumeration and end-to-end manifest contracts (stdlib only)."""
import cmath
import itertools
import math
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys

HPHI = str(Path(sys.argv[1]).resolve())
ROOT = Path("symmetry_sector_manifest")
ROOT.mkdir(exist_ok=True)
MPI = shlex.split(os.environ.get("MPIRUN", ""))
MASK = (1 << 64) - 1


def fnv(data):
    h = 14695981039346656037
    for byte in data:
        h = ((h ^ byte) * 1099511628211) & MASK
    return h


def manifest(path):
    rows = path.read_text().splitlines()
    values = dict(row.split("=", 1) for row in rows)
    assert len(rows) == len(values), "duplicate manifest key"
    assert values["format"] == "HPhiSymmetrySector version=1"
    return values


def run(path, label, args, layout="distributed", mpi=False, fail=False):
    env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)
    command = (MPI if mpi else []) + [HPHI] + args
    log = path / (label + ".log")
    with log.open("w") as stream:
        result = subprocess.run(command, cwd=str(path), env=env,
                                stdout=stream, stderr=subprocess.STDOUT, timeout=90)
    if (result.returncode != 0) != fail:
        raise AssertionError("{}: {}\n{}".format(label, command, log.read_text()))
    return log.read_text()


def write_def(path, name, count, rows):
    (path / name).write_text("====\nNItems {}\n====\n====\n====\n{}\n".format(count, "\n".join(rows)))


def write_group(path, length, k, reverse=False, metadata=True):
    order = list(range(length))
    if reverse:
        order.reverse()
    text = "# MomentumIndex {}\n".format(k) if metadata else ""
    text += "====\nNQPTrans {}\n====\n====\n====\n".format(length)
    for new, old in enumerate(order):
        ch = cmath.exp(-2j * math.pi * k * old / length)
        text += "{} {:.17g} {:.17g}\n".format(new, ch.real, ch.imag)
    for new, old in enumerate(order):
        for site in range(length):
            text += "{} {} {} 1\n".format(new, site, (site + old) % length)
    (path / "qptransidx.def").write_text(text)


def reference_sector(length, k, model):
    # Enumerate Fock configurations directly, preserving fermion permutation signs.
    filling = length // 2 if model != "Hubbard" else 1
    combinations = list(itertools.combinations(range(length), filling))
    if model == "Hubbard":
        configs = [tuple(sorted([2*i for i in up] + [2*i+1 for i in down]))
                   for up in combinations for down in combinations]
    else:
        configs = combinations
    seen = set()
    entries = []
    for occupied in configs:
        state = sum(1 << i for i in occupied)
        if state in seen:
            continue
        images = []
        for shift in range(length):
            moved = [2*((i//2 + shift) % length) + i % 2 if model == "Hubbard"
                     else (i+shift) % length for i in occupied]
            inversions = sum(a > b for pos, a in enumerate(moved) for b in moved[pos+1:])
            sign = 1 if model == "Spin" else (-1)**inversions
            images.append((sum(1 << i for i in moved), sign))
        orbit = {image for image, sign in images}
        seen.update(orbit)
        projector = sum(cmath.exp(-2j*math.pi*k*g/length)*sign
                        for g, (image, sign) in enumerate(images) if image == state)
        if abs(projector) > 1e-8:
            entries.append((min(orbit), len(orbit), length // len(orbit)))
    hashes = [fnv(b"hphi-sector-multiset-v1\0" + struct.pack("<QII", *entry))
              for entry in entries]
    xor = 0
    for h in hashes:
        xor ^= h
    return len(configs), len(entries), "hphi-sector-multiset-v1:{}:{:016x}:{:016x}".format(
        len(entries), xor, sum(hashes) & MASK)


def reference_group(length, k):
    data = b"hphi-group-fnv1a64-v1\0" + struct.pack("<II", length, length)
    for shift in range(length):  # For a chain, lexicographic order is shift order.
        for site in range(length):
            data += struct.pack("<Ii", (site+shift) % length, 1)
        ch = cmath.exp(-2j*math.pi*k*shift/length)
        data += struct.pack("<qq", round(1e10*ch.real), round(1e10*ch.imag))
    return "hphi-group-fnv1a64-v1:{:016x}".format(fnv(data))


def prepare(name, model, length):
    path = ROOT / name
    if path.exists():
        shutil.rmtree(str(path))
    path.mkdir()
    calc_model = {"Spin": 1, "SpinlessFermion": 7, "Hubbard": 0}[model]
    (path / "calcmod.def").write_text("CalcType 3\nCalcModel {}\nOutputMode 0\nOutputDataHead 1\n".format(calc_model))
    quantum = "2Sz 0\n" if model == "Spin" else (
        "Ncond {}\n".format(length//2) if model == "SpinlessFermion" else "Ncond 2\n2Sz 0\n")
    (path / "modpara.def").write_text(
        "====\nModel_Parameters 0\n====\n====\n====\nCDataFileHead trial\nCParaFileHead zqp\n====\n"
        "Nsite {}\n{}Lanczos_max 80\ninitial_iv -1\nexct 1\nLanczosEps 10\nLargeValue 50\n".format(length, quantum))
    write_def(path, "locspn.def", length if model == "Spin" else 0,
              ["{} {}".format(i, int(model == "Spin")) for i in range(length)])
    interaction = "Exchange exchange.def\n" if model == "Spin" else "Trans trans.def\n"
    if model == "Spin":
        write_def(path, "exchange.def", length,
                  ["{} {} 1.0".format(i, (i+1) % length) for i in range(length)])
    else:
        terms = []
        for i in range(length):
            j = (i+1) % length
            for spin in range(2 if model == "Hubbard" else 1):
                # Trans requires both real and imaginary coefficient columns.
                terms.extend(["{} {} {} {} -1.0 0.0".format(i, spin, j, spin),
                              "{} {} {} {} -1.0 0.0".format(j, spin, i, spin)])
        write_def(path, "trans.def", len(terms), terms)
        if model == "Hubbard":
            write_def(path, "coulombintra.def", length, ["{} 0.5".format(i) for i in range(length)])
            interaction += "CoulombIntra coulombintra.def\n"
    (path / "namelist.def").write_text("CalcMod calcmod.def\nModPara modpara.def\nLocSpin locspn.def\n" + interaction + "TransSym qptransidx.def\n")
    write_group(path, length, 1)
    return path


for model, length in [("Spin", 8), ("SpinlessFermion", 6), ("Hubbard", 4)]:
    path = prepare(model, model, length)
    full_dim, sector_dim, digest = reference_sector(length, 1, model)
    run(path, "serial", ["-e", "namelist.def"])
    file = path / "output/trial_symmetry_sector.dat"
    base = manifest(file)
    assert base["model"] == model and base["calc_type"] == "CG"
    assert int(base["full_dim"]) == full_dim and int(base["sector_dim"]) == sector_dim
    assert base["sector_digest"] == digest
    assert base["group_digest"] == reference_group(length, 1)
    if model == "Spin":
        assert base["fixed_2sz"] == "0" and "fixed_ne" not in base
    else:
        assert int(base["fixed_ne"]) == (2 if model == "Hubbard" else length//2)
        if model == "Hubbard":
            assert base["fixed_nup"] == base["fixed_ndown"] == "1"
    assert base["momentum_index"] == "1" and int(base["group_order"]) == length
    identity_keys = ["group_digest", "sector_digest", "hamiltonian_digest"]
    for layout in ["replicated", "distributed"]:
        run(path, layout, ["-e", "namelist.def"], layout, mpi=True)
        parallel = manifest(file)
        assert parallel["basis_layout"] == layout
        for key in identity_keys:
            assert parallel[key] == base[key], (
                model, layout, key, "serial", base[key], "parallel", parallel[key])
    write_group(path, length, 1, reverse=True, metadata=False)
    run(path, "renumbered", ["-e", "namelist.def"])
    reordered = manifest(file)
    assert "momentum_index" not in reordered
    assert all(reordered[key] == base[key] for key in identity_keys)
    write_group(path, length, length-1)
    run(path, "conjugate", ["-e", "namelist.def"])
    conjugate = manifest(file)
    assert conjugate["group_digest"] != base["group_digest"]
    assert conjugate["sector_digest"] == base["sector_digest"]
    assert conjugate["hamiltonian_digest"] == base["hamiltonian_digest"]
    write_group(path, length, 1)
    parameter = path / ("exchange.def" if model == "Spin" else "trans.def")
    parameter.write_text(parameter.read_text().replace("1.0", "0.7"))
    run(path, "new_coupling", ["-e", "namelist.def"])
    changed = manifest(file)
    assert changed["hamiltonian_digest"] != base["hamiltonian_digest"]
    assert changed["sector_digest"] == base["sector_digest"]
    assert changed["group_digest"] == base["group_digest"]
    if model == "Hubbard":
        interaction = path / "coulombintra.def"
        interaction.write_text(interaction.read_text().replace("0.5", "1.5"))
        run(path, "new_diagonal", ["-e", "namelist.def"])
        diagonal = manifest(file)
        assert diagonal["hamiltonian_digest"] != changed["hamiltonian_digest"]
        assert diagonal["sector_digest"] == base["sector_digest"]
    calcmod = path / "calcmod.def"
    calcmod.write_text(calcmod.read_text().replace("OutputDataHead 1", "OutputDataHead 0"))
    run(path, "no_prefix", ["-e", "namelist.def"])
    plain = path / "output/symmetry_sector.dat"
    assert manifest(plain)["sector_digest"] == base["sector_digest"]
    # A directory at the target path fails even when CI runs as root.
    plain.unlink()
    plain.mkdir()
    log = run(path, "write_failure", ["-e", "namelist.def"], mpi=True, fail=True)
    assert "could not write TransSym sector manifest" in log

# Standard-mode metadata propagation and dry-run exclusion.
path = ROOT / "standard"
path.mkdir(exist_ok=True)
if (path / "output").exists():
    shutil.rmtree(str(path / "output"))
(path / "stan.in").write_text("L=8\nmodel=Spin\nmethod=CG\nlattice=chain\noutputmode=none\n"
                             "Jx=1\nJy=1\nJz=0\n2Sz=0\nMomentumIndex=1\nexct=1\n")
run(path, "dry", ["-sdry", "stan.in"])
assert not list(path.glob("output/*symmetry_sector.dat"))
run(path, "standard", ["-e", "namelist.def"])
files = list(path.glob("output/*symmetry_sector.dat"))
assert len(files) == 1 and manifest(files[0])["momentum_index"] == "1"
# Ordinary runs do not acquire a manifest.
raw = ROOT / "raw"
raw.mkdir(exist_ok=True)
if (raw / "output").exists():
    shutil.rmtree(str(raw / "output"))
(raw / "stan.in").write_text((path / "stan.in").read_text().replace("MomentumIndex=1\n", ""))
run(raw, "raw", ["-s", "stan.in"], layout="replicated")
assert not list(raw.glob("output/*symmetry_sector.dat"))
print("sector manifest: independent Spin/fermion orbit digests, layouts, conjugates, coupling changes and I/O PASS")

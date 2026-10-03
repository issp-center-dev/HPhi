"""Raw TE I/O, first-step and restart contracts; no hidden MPI on negatives."""
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import sys

import numpy as np

HPHI = str(Path(sys.argv[1]).resolve())
MPI = shlex.split(os.environ.get("MPIRUN", ""))
ROOT = Path("te_contract")
if ROOT.exists():
    shutil.rmtree(ROOT)
ROOT.mkdir()


def run(path, label, failure=None, mpi=True):
    with (path / (label + ".log")).open("w") as fp:
        result = subprocess.run((MPI if mpi else []) + [HPHI, "-e", "namelist.def"],
                                cwd=str(path), stdout=fp, stderr=subprocess.STDOUT, timeout=60)
    text = (path / (label + ".log")).read_text()
    if failure:
        assert result.returncode != 0 and failure in text, text
    else:
        assert result.returncode == 0, text


def check_final_extension(source, family=None):
    """A final file from six rows must resume at row six on a longer grid."""
    path = source / "final_extension"
    path.mkdir()
    (path / "output").mkdir()
    for p in source.glob("*.def"):
        shutil.copyfile(p, path / p.name)
    for name in seed_names:
        shutil.copyfile(source / "output" / name, path / "output" / name)
    for rank, seed in enumerate(seeds):
        final = source / "output" / seed.name
        assert struct.unpack("=i", final.read_bytes()[:4])[0] == 6
        shutil.copyfile(final, path / "output/extend_rank_{}.dat".format(rank))
    p = path / "modpara.def"
    p.write_text(re.sub(r"^\s*Lanczos_max\s+\d+", "Lanczos_max 9", p.read_text(), flags=re.M))
    rows = "====\nNTimeSteps 9\n====\n====\n====\n"
    for step in range(9):
        rows += "{} {}\n".format(.025*step, int(family is not None))
        if family:
            indices = "0 0 0 0" if family == "TEOneBody" else "0 0 0 0 1 0 1 0"
            rows += indices + " {} 0\n".format(.2*step)
    (path / "times.def").write_text(rows)
    run(path, "uninterrupted")
    expected = [(path / "output" / seed.name).read_bytes() for seed in seeds]
    observables = {key: np.loadtxt(path / "output" / (key + ".dat"), ndmin=2)
                   for key in ("SS", "Norm", "Flct")}
    shutil.copytree(path / "output", path / "uninterrupted_output")
    p = path / "namelist.def"
    p.write_text(p.read_text().replace("seed%literal_eigenvec_0", "extend"))
    p = path / "calcmod.def"
    p.write_text(re.sub(r"^ReStart\s+\d+", "ReStart 3", p.read_text(), flags=re.M))
    run(path, "extended")
    for seed, saved in zip(seeds, expected):
        actual = (path / "output" / seed.name).read_bytes()
        assert struct.unpack("=i", actual[:4])[0] == 9
        assert actual[:12] == saved[:12]
        np.testing.assert_allclose(np.frombuffer(actual[12:], dtype=complex),
                                   np.frombuffer(saved[12:], dtype=complex), atol=3e-13, rtol=0)
    for key, reference in observables.items():
        actual = np.loadtxt(path / "output" / (key + ".dat"), ndmin=2)
        np.testing.assert_allclose(actual, reference[6:], atol=5e-12, rtol=0)


(ROOT / "stan.in").write_text('L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\nexct=1\noutputmode="None"\n')
with (ROOT / "generate.log").open("w") as fp:
    result = subprocess.run([HPHI, "-sdry", "stan.in"], cwd=str(ROOT), stdout=fp, stderr=subprocess.STDOUT)
assert result.returncode == 0, (ROOT / "generate.log").read_text()
calc = (ROOT / "calcmod.def").read_text()
calc = re.sub(r"^OutputEigenVec\s+\d+", "OutputEigenVec 1", calc, flags=re.M)
(ROOT / "calcmod.def").write_text(calc)
run(ROOT, "cg")
seeds = sorted((ROOT / "output").glob("zvo_eigenvec_0_rank_*.dat"),
               key=lambda p: int(p.stem.split("rank_")[1]))
assert seeds
for p in seeds:
    shutil.copyfile(p, ROOT / "output" / p.name.replace("zvo_", "seed%literal_"))
seed_names = [p.name.replace("zvo_", "seed%literal_") for p in seeds]
# Quench the transverse exchange while retaining the longitudinal interaction.
p = ROOT / "exchange.def"
lines = p.read_text().splitlines()
for i in range(5, len(lines)):
    row = lines[i].split()
    if row:
        row[-1] = str(float(row[-1]) * 1.3)
        lines[i] = " ".join(row)
p.write_text("\n".join(lines) + "\n")
calc = re.sub(r"^CalcType\s+\d+", "CalcType 4", calc, flags=re.M)
calc = re.sub(r"^InputEigenVec\s+\d+", "InputEigenVec 1", calc, flags=re.M)
(ROOT / "calcmod.def").write_text(calc)
mod = (ROOT / "modpara.def").read_text()
mod = re.sub(r"^\s*Lanczos_max\s+\d+", "Lanczos_max 6", mod, flags=re.M)
mod = "\n".join(line for line in mod.splitlines() if not line.strip().startswith(("ExpandCoef", "OutputInterval"))) + "\n"
(ROOT / "modpara.def").write_text(mod + "ExpandCoef 8\nOutputInterval 1\n")
names = (ROOT / "namelist.def").read_text()
names = "\n".join(line for line in names.splitlines() if not line.lstrip().startswith("SpectrumVec"))
names += "\nSpectrumVec seed%literal_eigenvec_0\nTEOneBody times.def\n"
(ROOT / "namelist.def").write_text(names)
(ROOT / "times.def").write_text("====\nNTimeSteps 6\n====\n====\n====\n" +
                                "".join("{} 0\n".format(.025*i) for i in range(6)))
run(ROOT, "uninterrupted")
expected = [p.read_bytes() for p in seeds]  # legacy final alias has index zero
observables = {family: np.loadtxt(ROOT / "output" / (family + ".dat")) for family in ("SS", "Norm", "Flct")}
shutil.copytree(ROOT / "output", ROOT / "uninterrupted_output")
check_final_extension(ROOT)
# A periodic file after step 2 must encode NEXT step 3, just like final files.
for rank in range(len(seeds)):
    source = ROOT / "output/zvo_eigenvec_2_rank_{}.dat".format(rank)
    assert struct.unpack("=i", source.read_bytes()[:4])[0] == 3
    shutil.copyfile(source, ROOT / "output/restart_rank_{}.dat".format(rank))
(ROOT / "namelist.def").write_text(names.replace("seed%literal_eigenvec_0", "restart"))
(ROOT / "calcmod.def").write_text(re.sub(r"^ReStart\s+\d+", "ReStart 3", calc, flags=re.M))
run(ROOT, "restart")
for rank, saved in enumerate(expected):
    actual = seeds[rank].read_bytes()
    assert actual[:12] == saved[:12]
    np.testing.assert_allclose(np.frombuffer(actual[12:], dtype=complex),
                               np.frombuffer(saved[12:], dtype=complex), atol=3e-13, rtol=0)
for family, reference in observables.items():
    actual = np.loadtxt(ROOT / "output" / (family + ".dat"), ndmin=2)
    np.testing.assert_allclose(actual, reference[3:], atol=5e-12, rtol=0)
shutil.copytree(ROOT / "output", ROOT / "restart_output")

# Single-rank damaged files must fail collectively before creating thermal output.
(ROOT / "namelist.def").write_text(names)
(ROOT / "calcmod.def").write_text(calc)
victim = ROOT / "output" / seed_names[-1]
saved = victim.read_bytes()
for label, data in [("short_header", saved[:2]), ("short_payload", saved[:-1]),
                    ("trailing", saved+b"x")]:
    victim.write_bytes(data)
    run(ROOT, label, "missing, truncated or incompatible TE Inputvector")
victim.unlink()
run(ROOT, "missing", "missing, truncated or incompatible TE Inputvector")
victim.write_bytes(saved)
restart_victim = ROOT / "output/restart_rank_{}.dat".format(len(seeds)-1)
restart_saved = restart_victim.read_bytes()
bad_step = bytearray(restart_saved)
bad_step[:4] = struct.pack("=i", 6)
restart_victim.write_bytes(bad_step)
(ROOT / "namelist.def").write_text(names.replace("seed%literal_eigenvec_0", "restart"))
(ROOT / "calcmod.def").write_text(re.sub(r"^ReStart\s+\d+", "ReStart 3", calc, flags=re.M))
run(ROOT, "out_of_range", "TE restart step must be consistent")
restart_victim.write_bytes(restart_saved)
(ROOT / "namelist.def").write_text(names)
(ROOT / "calcmod.def").write_text(calc)
(ROOT / "modpara.def").write_text(mod + "OutputInterval 1\n")
run(ROOT, "missing_order", "TE requires positive ExpandCoef")
(ROOT / "modpara.def").write_text(mod + "ExpandCoef 8\nOutputInterval 1\n")
blocked = ROOT / "output/zvo_eigenvec_1_rank_{}.dat".format(len(seeds)-1)
blocked.unlink()
blocked.mkdir()
run(ROOT, "write_failure", "failed to write TE vector")
blocked.rmdir()
for keyword, bad_values in [("Lanczos_max", ["0", "-1", "1.5", "nan", "1e100"]),
                            ("OutputInterval", ["0", "-1", "1.5", "nan", "1e100"])]:
    for bad in bad_values:
        badmod = re.sub(r"^\s*" + keyword + r"\s+[^\n]+", keyword + " " + bad,
                        mod + "ExpandCoef 8\nOutputInterval 1\n", flags=re.M)
        (ROOT / "modpara.def").write_text(badmod)
        # Definition-only guards do not need MPI.
        run(ROOT, "bad_{}_{}".format(keyword, bad), keyword + " must be a positive integer", mpi=False)
(ROOT / "modpara.def").write_text(mod + "ExpandCoef 8\nOutputInterval 1\n")

for family in ("TEOneBody", "TETwoBody"):
    path = ROOT / family
    path.mkdir()
    (path / "output").mkdir()
    for p in ROOT.glob("*.def"):
        shutil.copyfile(p, path / p.name)
    for name in seed_names:
        shutil.copyfile(ROOT / "output" / name, path / "output" / name)
    (path / "calcmod.def").write_text(calc)
    dynamic_names = names.replace("TEOneBody times.def", family + " times.def")
    (path / "namelist.def").write_text(dynamic_names)
    rows = "====\nNTimeSteps 6\n====\n====\n====\n"
    for step in range(6):
        rows += "{} 1\n".format(.025*step)
        indices = "0 0 0 0" if family == "TEOneBody" else "0 0 0 0 1 0 1 0"
        rows += indices + " {} 0\n".format(.2*step)
    (path / "times.def").write_text(rows)
    run(path, "uninterrupted")
    check_final_extension(path, family)
    reference = [(path / "output" / p.name).read_bytes() for p in seeds]
    ss = np.loadtxt(path / "output/SS.dat")
    for rank in range(len(seeds)):
        shutil.copyfile(path / "output/zvo_eigenvec_2_rank_{}.dat".format(rank),
                        path / "output/restart_rank_{}.dat".format(rank))
    shutil.copytree(path / "output", path / "uninterrupted_output")
    (path / "namelist.def").write_text(dynamic_names.replace("seed%literal_eigenvec_0", "restart"))
    (path / "calcmod.def").write_text(re.sub(r"^ReStart\s+\d+", "ReStart 3", calc, flags=re.M))
    run(path, "restart")
    for rank, saved in enumerate(reference):
        actual = (path / "output" / seeds[rank].name).read_bytes()
        np.testing.assert_allclose(np.frombuffer(actual[12:], dtype=complex),
                                   np.frombuffer(saved[12:], dtype=complex), atol=3e-13, rtol=0)
    np.testing.assert_allclose(np.loadtxt(path / "output/SS.dat"), ss[3:], atol=5e-12, rtol=0)
print("Raw TE: literal filenames, collective I/O guards, periodic restart and final extension PASS")

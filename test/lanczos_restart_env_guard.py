"""Ground-state Lanczos restart validates its MPI/OpenMP environment."""
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


HPHI = str(Path(sys.argv[1]).resolve())
MPI = shlex.split(os.environ.get("MPIRUN", ""))
ROOT = Path("lanczos_restart_env_guard")
if ROOT.exists():
    shutil.rmtree(ROOT)
ROOT.mkdir()
STAN = (
    'L=8\nmodel="Spin"\nmethod="Lanczos"\nlattice="chain"\nJ=1\n2Sz=0\n'
    'Lanczos_max={}\nRestart="{}"\noutputmode="None"\n'
)


def execute(label, lanczos_max, restart, threads):
    (ROOT / "stan.in").write_text(STAN.format(lanczos_max, restart))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    log = ROOT / (label + ".log")
    with log.open("w") as fp:
        result = subprocess.run(
            MPI + [HPHI, "-s", "stan.in"],
            cwd=str(ROOT),
            env=env,
            stdout=fp,
            stderr=subprocess.STDOUT,
            timeout=300,
        )
    return result.returncode, log.read_text()


def run(label, lanczos_max, restart, threads):
    status, text = execute(label, lanczos_max, restart, threads)
    assert status == 0, text
    return text


def strip_meta():
    tm = ROOT / "output/zvo_TMComponents.dat"
    tm.write_text(
        "\n".join(
            line
            for line in tm.read_text().splitlines()
            if not line.lstrip().startswith("#")
        )
        + "\n"
    )


run("out", 10, "restart_out", 1)
tm_text = (ROOT / "output/zvo_TMComponents.dat").read_text()
match = re.search(r"^# nproc (\d+) nthreads 1$", tm_text, flags=re.M)
assert match, tm_text
saved_nproc = int(match.group(1))
assert saved_nproc == max(1, len(MPI) and int(MPI[-1]) or 1), (saved_nproc, MPI)
run("same", 1000, "restart_in", 1)

run("out", 10, "restart_out", 1)
status, mismatch_text = execute("mismatch", 1000, "restart_in", 2)
actual_threads = re.search(r"OpenMP threads\s*:\s*(\d+)", mismatch_text)
if actual_threads and int(actual_threads.group(1)) != 1:
    assert status != 0 and "TMComponents were written with nthreads=1" in mismatch_text, mismatch_text
else:
    # Builds without OpenMP always report one thread; the metadata path is still
    # covered here, while OpenMP-enabled CI exercises the rejection above.
    assert status == 0, mismatch_text

run("out", 10, "restart_out", 1)
strip_meta()
legacy_text = run("legacy", 1000, "restart_in", 1)
assert "Warning: TMComponents has no environment line" in legacy_text, legacy_text
print("ok")

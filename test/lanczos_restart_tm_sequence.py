"""Restarted Lanczos reproduces the uninterrupted tridiagonal sequence."""
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


HPHI = str(Path(sys.argv[1]).resolve())
MPI = shlex.split(os.environ.get("MPIRUN", ""))
ROOT = Path("lanczos_restart_tm_sequence")
shutil.rmtree(ROOT, ignore_errors=True)
ROOT.mkdir()
STAN = (
    'L=10\nmodel="Spin"\nmethod="Lanczos"\nlattice="chain"\nJ=1\n2Sz=0\n'
    'Lanczos_max={}\nLanczosEps=16\nRestart="{}"\noutputmode="None"\n'
)


def run(path, label, lanczos_max, restart):
    path.mkdir(exist_ok=True)
    (path / "stan.in").write_text(STAN.format(lanczos_max, restart))
    log = path / (label + ".log")
    with log.open("w") as fp:
        result = subprocess.run(
            MPI + [HPHI, "-s", "stan.in"],
            cwd=str(path),
            stdout=fp,
            stderr=subprocess.STDOUT,
            timeout=300,
        )
    assert result.returncode == 0, log.read_text()


def tm(path):
    lines = [
        line
        for line in (path / "output/zvo_TMComponents.dat").read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return int(lines[0]), np.array(
        [[float(value) for value in line.split()] for line in lines[2:]]
    )


reference = ROOT / "ref"
run(reference, "ref", 20, "restart_out")
reference_steps, reference_components = tm(reference)

continued = ROOT / "cont"
run(continued, "first", 10, "restart_out")
run(continued, "second", 10, "restart")
continued_steps, continued_components = tm(continued)

assert continued_steps == reference_steps == 20, (continued_steps, reference_steps)
assert continued_components.shape == reference_components.shape, (
    continued_components.shape,
    reference_components.shape,
)
np.testing.assert_allclose(
    continued_components, reference_components, rtol=0, atol=2e-9
)
print("ok")

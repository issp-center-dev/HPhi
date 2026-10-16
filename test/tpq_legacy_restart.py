"""Legacy tmpvec continuation restores cTPQ inverse temperature and resets firstnorm per set."""
import os
import shlex
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np


HPHI = str(Path(sys.argv[1]).resolve())
MPI = shlex.split(os.environ.get("MPIRUN", ""))
ROOT = Path("tpq_legacy_restart")
if ROOT.exists():
    shutil.rmtree(ROOT)
ROOT.mkdir()


def write_stan(path, method, lanczos_max, restart, numave=2):
    path.mkdir(exist_ok=True)
    (path / "stan.in").write_text(
        'L=8\nmodel="Spin"\nmethod="{}"\nlattice="chain"\nJ=1\n2Sz=0\nNumAve={}\n'
        'Lanczos_max={}\nLargeValue=4\nExpecInterval=1000\nRestart="{}"\noutputmode="None"\n'
        .format(method, numave, lanczos_max, restart)
    )


def run(path, label):
    with (path / (label + ".log")).open("w") as fp:
        result = subprocess.run(
            MPI + [HPHI, "-s", "stan.in"],
            cwd=str(path),
            stdout=fp,
            stderr=subprocess.STDOUT,
            timeout=300,
        )
    assert result.returncode == 0, (path / (label + ".log")).read_text()


def rows(path, name):
    return np.loadtxt(path / "output" / name, ndmin=2)


def make_legacy_fixture(path):
    """No-op with the PR-1 writer; explicitly convert all manifest rows with the PR-2 writer."""
    output = path / "output"
    manifest = output / "zvo_restart_manifest.dat"
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            fields = line.split()
            if not fields or fields[0] != "set":
                continue
            set_index, ranks, generation = int(fields[1]), int(fields[3]), int(fields[5])
            for rank in range(ranks):
                data = (
                    output
                    / "zvo_restart_set{}_comp0_gen{}_rank_{}.dat".format(
                        set_index, generation, rank
                    )
                ).read_bytes()
                words = struct.unpack_from("<40Q", data)
                dimension, step = words[15], words[24]
                assert words[1] == 3 and len(data) == 320 + 16 * dimension
                values = np.frombuffer(data, dtype="<c16", offset=320).astype(np.complex128)
                (output / "tmpvec_set{}_rank_{}.dat".format(set_index, rank)).write_bytes(
                    struct.pack("=iQ", step, dimension) + bytes(16) + values.tobytes()
                )
        manifest.unlink()
        for checkpoint in output.glob("zvo_restart_set*"):
            checkpoint.unlink()
    assert list(output.glob("tmpvec_set*_rank_*.dat")), output


for method in ("cTPQ", "TPQ"):
    reference = ROOT / (method + "_ref")
    write_stan(reference, method, 10, "none")
    run(reference, "ref")

    continued = ROOT / (method + "_cont")
    write_stan(continued, method, 6, "restart_out")
    run(continued, "first")
    make_legacy_fixture(continued)
    write_stan(continued, method, 10, "restart_in")
    run(continued, "second")
    for set_index in range(2):
        expected = rows(reference, "SS_rand{}.dat".format(set_index))
        actual = rows(continued, "SS_rand{}.dat".format(set_index))
        np.testing.assert_allclose(
            actual[actual[:, -1] >= 6][:, :2],
            expected[expected[:, -1] >= 6][:, :2],
            rtol=0,
            atol=1e-10,
        )
        if method == "cTPQ":
            np.testing.assert_allclose(
                actual[actual[:, -1] >= 6][:, 0],
                actual[actual[:, -1] >= 6][:, -1] / 4.0,
                rtol=0,
                atol=1e-12,
            )

    mixed = ROOT / (method + "_mixed")
    write_stan(mixed, method, 6, "restart_out")
    run(mixed, "first")
    make_legacy_fixture(mixed)
    for checkpoint in (mixed / "output").glob("tmpvec_set0_rank_*.dat"):
        checkpoint.unlink()
    write_stan(mixed, method, 10, "restart")
    run(mixed, "second")
    norm1 = rows(mixed, "Norm_rand1.dat")
    assert np.all(norm1[norm1[:, -1] >= 6][:, 2] == 0.0), norm1
    assert "Start to calculate in normal procedure" in (mixed / "second.log").read_text()

print("ok")

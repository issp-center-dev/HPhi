"""Shared helpers for restart tests: input generation, runs, manifest and v3 header access."""
import os, re, shlex, shutil, struct, subprocess, sys
from pathlib import Path
import numpy as np

MPI = shlex.split(os.environ.get("MPIRUN", ""))
V3_WORDS = 40
MAGIC = 0x0a31565349485048

def mpi_size():
    """Return the process count requested by MPIRUN, or one when unset."""
    for i, arg in enumerate(MPI[:-1]):
        if arg in ("-np", "-n"):
            return int(MPI[i + 1])
    return 1

def launcher(ranks):
    """Return MPIRUN with its process count replaced by *ranks*.

    Rank-redistribution tests execute both the writer and reader from one
    CTest process.  A serial run needs no launcher; a non-serial run without
    an explicit ``-np``/``-n`` launcher cannot be expressed portably and is
    skipped with CTest's conventional return code.
    """
    if ranks < 1:
        raise ValueError("MPI rank count must be positive")
    if not MPI:
        if ranks == 1:
            return []
        print("Skipping: MPIRUN is required for {} ranks".format(ranks), flush=True)
        sys.exit(77)
    command = list(MPI)
    for i, arg in enumerate(command[:-1]):
        if arg in ("-np", "-n"):
            command[i + 1] = str(ranks)
            return command
    print("Skipping: could not replace -np/-n in MPIRUN={!r}".format(MPI), flush=True)
    sys.exit(77)

def hphi():
    return str(Path(sys.argv[1]).resolve())

def probe():
    return str(Path(hphi()).with_name("unittest_symmetry_sector_probe"))

def fresh_dir(name):
    root = Path(name)
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    return root

def generate(path, stan_text):
    """Expert-mode files from a Standard-mode input via -sdry (serial)."""
    path.mkdir(parents=True, exist_ok=True)
    (path / "stan.in").write_text(stan_text)
    with (path / "generate.log").open("w") as fp:
        r = subprocess.run([hphi(), "-sdry", "stan.in"], cwd=str(path), stdout=fp, stderr=subprocess.STDOUT)
    assert r.returncode == 0, (path / "generate.log").read_text()

def set_key(path, filename, key, value):
    """Set or append `key value` in a def file (calcmod.def / modpara.def)."""
    p = path / filename
    text = p.read_text()
    new, n = re.subn(r"^(\s*)" + re.escape(key) + r"\s+\S+.*$", r"\g<1>" + key + " " + str(value), text, flags=re.M)
    p.write_text(new if n else text.rstrip("\n") + "\n" + key + " " + str(value) + "\n")

def drop_key(path, filename, key):
    p = path / filename
    p.write_text("\n".join(l for l in p.read_text().splitlines() if not l.strip().startswith(key)) + "\n")

def run(path, label, binary=None, env=None, fail=None, mpi=True, timeout=600,
        namelist="namelist.def"):
    cmd = (MPI if mpi else []) + [binary or hphi(), "-e", namelist]
    e = dict(os.environ); e.update(env or {})
    log = path / (label + ".log")
    with log.open("w") as fp:
        r = subprocess.run(cmd, cwd=str(path), env=e, stdout=fp, stderr=subprocess.STDOUT, timeout=timeout)
    text = log.read_text()
    if fail is not None:
        assert r.returncode != 0 and fail in text, text
    else:
        assert r.returncode == 0, text
    return text

def manifest(path, head="zvo"):
    p = path / "output" / (head + "_restart_manifest.dat")
    if not p.exists():
        return None
    rows = {}
    meta = {}
    for line in p.read_text().splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "set":
            rows[int(t[1])] = dict(ranks=int(t[3]), gen=int(t[5]), step=int(t[7]), ncomp=int(t[9]))
        else:
            meta[t[0]] = t[1]
    return meta, rows

def v3_header(path):
    data = Path(path).read_bytes()
    words = struct.unpack("<{}Q".format(V3_WORDS), data[:8 * V3_WORDS])
    assert words[0] == MAGIC and words[1] == 3, words[:2]
    payload = np.frombuffer(data[8 * V3_WORDS:], dtype="<c16")
    assert len(payload) == words[15], (len(payload), words[15])
    return words, payload

def restart_files(path, head, k, comp, gen):
    return sorted((path / "output").glob("{}_restart_set{}_comp{}_gen{}_rank_*.dat".format(head, k, comp, gen)),
                  key=lambda p: int(p.stem.split("rank_")[1]))

def v3_to_legacy_tmpvec(path, head, set_index):
    """Convert one v3 TPQ set into same-rank legacy tmpvec files."""
    found = manifest(path, head)
    assert found is not None and set_index in found[1], (path, head, set_index)
    row = found[1][set_index]
    files = restart_files(path, head, set_index, 0, row["gen"])
    assert len(files) == row["ranks"], (files, row)
    for rank, source in enumerate(files):
        words, payload = v3_header(source)
        assert words[23] == set_index and words[24] == row["step"]
        target = path / "output" / "tmpvec_set{}_rank_{}.dat".format(set_index, rank)
        target.write_bytes(struct.pack("=iQ", row["step"], words[15]) +
                           bytes(16) + payload.tobytes())

def rows(path, name):
    return np.loadtxt(path / "output" / name, ndmin=2)

def sorted_by_key(table, key_columns):
    order = np.lexsort([table[:, c] for c in reversed(key_columns)])
    return table[order]

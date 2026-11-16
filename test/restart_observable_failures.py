"""Actual output-open and checked-close failures at solver observable boundaries."""
import re
import shutil

import restart_common as rc

# One family at a time ensures each failing open reaches its own entry point.
FAMILIES = {
    "one": ("OneBodyG", "cisajs", "0 0 0 0"),
    "two": ("TwoBodyG", "cisajscktalt", "0 0 0 0 1 0 1 0"),
    "three": ("ThreeBodyG", "ThreeBody", "0 0 0 0 " * 3),
    "four": ("FourBodyG", "FourBody", "0 0 0 0 " * 4),
    "six": ("SixBodyG", "SixBody", "0 0 0 0 " * 6),
    "nbody": ("NBodyG", "NBodyG", "2 0 0 0 0 1 0 1 0"),
    "anomalous": ("AnomalousG", "AnomalousG", "0 0 0 0 1"),
}


def snapshot(path):
    result = {p.name: p.read_bytes() for p in
              (path / "output").glob("zvo_restart*")}
    assert result
    return result


def agreed(text, prefix="Green output agreed failure rank="):
    assert "MPI_ERR" not in text, text
    ranks = {int(r) for r in re.findall(re.escape(prefix) + r"([0-9]+)", text)}
    assert ranks == set(range(rc.mpi_size())), (ranks, rc.mpi_size(), text)


def aggregate_inner_open_failures(base):
    """Fail the inner aggregate opener after recovery skipped initialization."""
    for rank in sorted({0, rc.mpi_size() - 1}):
        bad = base.with_name(base.name + "_inner_open_rank{}".format(rank))
        shutil.copytree(base, bad)
        before = snapshot(bad)
        text = rc.run(bad, "inner_open_status", binary=rc.probe(), timeout=60,
                      fail="correlation output open failed",
                      env={"HPHI_TEST_FD_FAIL": "green_open",
                           "HPHI_TEST_FD_FAIL_RANK": str(rank)})
        agreed(text)
        assert snapshot(bad) == before


def exercise(seed, family, aggregate, next_step, te=False):
    keyword, stem, row = FAMILIES[family]
    name = "{}_{}_{}".format(seed.name, family, "aggregate" if aggregate else "split")
    base = seed.with_name(name)
    shutil.copytree(seed, base)
    rc.set_key(base, "calcmod.def", "ReStart", 2)
    if te:
        rc.set_key(base, "calcmod.def", "InputEigenVec", 0)
    rc.set_key(base, "calcmod.def", "OutputGreenFormat", int(aggregate))
    rc.set_key(base, "modpara.def", "ExpecInterval", 1)
    rc.set_key(base, "modpara.def", "Lanczos_max", next_step + 2)
    rc.set_key(base, "modpara.def", "CheckpointInterval", 0 if te else 1)
    names = base / "namelist.def"
    keys = tuple(v[0] for v in FAMILIES.values())
    names.write_text("\n".join(line for line in names.read_text().splitlines()
                               if not line.startswith(keys)) +
                     "\n{} requested.def\n".format(keyword))
    (base / "requested.def").write_text("====\nNRequested 1\n====\n====\n====\n" + row + "\n")
    suffix = "te" if te else "tpq"
    filename = ("zvo_{}_{}.dat".format(stem, suffix) if aggregate else
                "zvo_{}_{}step{}.dat".format(stem, "" if te else "set0", next_step))
    # Keep older successful observations: recovery sees a missing new file,
    # then the real libc fopen encounters a directory at the write boundary.
    # An aggregate filename is also read by recovery, so omit its recovery
    # entry by starting a fresh ReStart=1 while retaining the saved checkpoint.
    if aggregate:
        rc.set_key(base, "calcmod.def", "ReStart", 1)
        if te:
            rc.set_key(base, "calcmod.def", "InputEigenVec", 1)
    for binary_label, binary in (("production", rc.hphi()), ("probe", rc.probe())):
        bad = base.with_name(name + "_open_" + binary_label)
        shutil.copytree(base, bad)
        target = bad / "output" / filename
        if target.exists():
            target.rename(target.with_name(target.name + ".retained"))
        target.mkdir()
        before = snapshot(bad)
        text = rc.run(bad, "actual_open", binary=binary, timeout=60,
                      fail="correlation output open failed")
        assert "MPI_ERR" not in text, text
        if binary_label == "probe":
            agreed(text)
        assert snapshot(bad) == before
    # Close errors are injected after actual observable evaluation, on root
    # and the last rank. They test propagation, not libc fault interposition.
    # Resume avoids aggregate initialization so this reaches the helper close.
    rc.set_key(base, "calcmod.def", "ReStart", 2)
    if te:
        rc.set_key(base, "calcmod.def", "InputEigenVec", 0)
    if aggregate:
        # Distinct from the actual libc aggregate initialization failure above.
        aggregate_inner_open_failures(base)
    for rank in sorted({0, rc.mpi_size() - 1}):
        bad = base.with_name(name + "_close_rank{}".format(rank))
        shutil.copytree(base, bad)
        before = snapshot(bad)
        text = rc.run(bad, "close_status", binary=rc.probe(), timeout=60,
                      fail="Green output agreed failure" if te else "driver operation failed",
                      env={"HPHI_TEST_FD_FAIL": "green_close",
                           "HPHI_TEST_FD_FAIL_RANK": str(rank),
                           # TE closes SS/Norm/Flct before its first Green.
                           "HPHI_TEST_FD_FAIL_SKIP": "3" if te else "0"})
        agreed(text)
        assert snapshot(bad) == before


def tpq_failures(root, method):
    for model in ("SpinGC", "HubbardGC"):
        seed = root / ("output_failure_" + model)
        stan = ('L={}\nmodel="{}"\nmethod="{}"\nlattice="chain"\n'
                'NumAve=1\nLanczos_max=3\nLargeValue=4\n'
                'outputmode="None"\n').format(8 if model == "SpinGC" else 4, model, method)
        stan += "J=1\n" if model == "SpinGC" else "t=1\nU=4\n"
        rc.generate(seed, stan)
        rc.set_key(seed, "calcmod.def", "ReStart", 1)
        rc.run(seed, "saved")
        for family in (list(FAMILIES)[:-1] if model == "SpinGC" else ["anomalous"]):
            for aggregate in (False, True):
                exercise(seed, family, aggregate, 3)
    # Real root fopen failures for all physical-data streams, including their
    # aggregate initialization (which precedes the first numerical operation).
    for aggregate in (False, True):
        for family in ("SS", "Norm", "Flct"):
            bad = root / "physical_{}_{}".format(family, int(aggregate))
            shutil.copytree(seed, bad)
            rc.set_key(bad, "calcmod.def", "OutputGreenFormat", int(aggregate))
            filename = family + ("_tpq.dat" if aggregate else "_rand0.dat")
            target = bad / "output" / filename
            if target.exists():
                target.rename(target.with_name(target.name + ".retained"))
            target.mkdir()
            before = snapshot(bad)
            text = rc.run(bad, "actual_physical_open", binary=rc.probe(), timeout=60,
                          fail="output open failed")
            if aggregate:
                agreed(text)
            else:
                solver = "cTPQ" if method == "cTPQ" else "mTPQ"
                agreed(text, solver + " agreed failure stage=output_open rank=")
            assert snapshot(bad) == before

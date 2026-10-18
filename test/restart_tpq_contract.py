"""mTPQ/cTPQ continuation: per-set and shared outputs match a full run."""
import os
import re
import shutil
import struct
import sys

import numpy as np

import restart_common as rc

MODE = sys.argv[2] if len(sys.argv) > 2 else "tpq"
METHOD = {"tpq": "TPQ", "ctpq": "cTPQ"}[MODE]
TINY_STEP = {"tpq": 2, "ctpq": 1}[MODE]
ROOT = rc.fresh_dir("restart_tpq_contract_" + MODE)
STAN = ('L=8\nmodel="Spin"\nmethod="{}"\nlattice="chain"\nJ=1\n2Sz=0\n'
        'NumAve=2\nLanczos_max=20\nLargeValue=4\nExpecInterval=1000\n'
        'outputmode="None"\n').format(METHOD)


def add_observations(path):
    (path / "greenone.def").write_text(
        "===============================\nNCisAjs 2\n===============================\n"
        "======== Green functions ======\n===============================\n"
        "0 0 0 0\n1 1 1 1\n")
    namelist = path / "namelist.def"
    namelist.write_text(namelist.read_text() + "OneBodyG greenone.def\n")
    # cTPQ reaches step 1 through its regular loop, so request every step to
    # cover that method-specific Green boundary. mTPQ's bootstrap step 1 does
    # not emit Green even when every step is requested.
    rc.set_key(path, "modpara.def", "ExpecInterval",
               1 if MODE == "ctpq" else 2)


def prepare(name, restart, numave=2, lanczos_max=20, interval=0,
            aggregate=False, observations=False):
    path = ROOT / name
    rc.generate(path, STAN)
    rc.set_key(path, "calcmod.def", "ReStart", restart)
    rc.set_key(path, "modpara.def", "NumAve", numave)
    rc.set_key(path, "modpara.def", "Lanczos_max", lanczos_max)
    if interval:
        rc.set_key(path, "modpara.def", "CheckpointInterval", interval)
    if aggregate:
        rc.set_key(path, "calcmod.def", "OutputGreenFormat", 1)
    if observations:
        add_observations(path)
    return path


def compare_set_files(a, b, sets):
    for set_index in sets:
        for name in ("SS_rand{}.dat", "Norm_rand{}.dat", "Flct_rand{}.dat"):
            np.testing.assert_allclose(rc.rows(b, name.format(set_index)),
                                       rc.rows(a, name.format(set_index)), rtol=0, atol=1e-10)


def compare_aggregates(a, b):
    for name in ("SS_tpq.dat", "Norm_tpq.dat", "Flct_tpq.dat"):
        x = rc.sorted_by_key(rc.rows(a, name), (0, 1))
        y = rc.sorted_by_key(rc.rows(b, name), (0, 1))
        np.testing.assert_allclose(y, x, rtol=0, atol=1e-10)


def assert_ctpq_beta_grid(path, aggregate=False, sets=(0, 1)):
    if MODE != "ctpq":
        return
    def data_name(name):
        return name if (path / "output" / name).exists() else "zvo_" + name
    if aggregate:
        table = rc.rows(path, data_name("SS_tpq.dat"))
        np.testing.assert_allclose(table[:, 2], table[:, 1] / 4.0,
                                   rtol=0, atol=1e-14)
    else:
        for set_index in sets:
            table = rc.rows(path, data_name("SS_rand{}.dat".format(set_index)))
            np.testing.assert_allclose(table[:, 0], table[:, -1] / 4.0,
                                       rtol=0, atol=1e-14)


def word_double(word):
    return struct.unpack("<d", struct.pack("<Q", word))[0]


def sorted_observations(path, name="zvo_cisajs_tpq.dat"):
    table = rc.rows(path, name)
    key = tuple(range(table.shape[1] - 2))
    return rc.sorted_by_key(table, key)


def is_physical_series(name):
    """SS/Norm/Flct and aggregate/split Green, with any data-head prefix."""
    data = re.search(r"(?:^|_)(?:SS|Norm|Flct)_(?:rand[0-9]+|tpq)[.]dat$", name)
    split_green = re.search(r"_set[0-9]+step[0-9]+[.]dat$", name)
    aggregate_green = name.endswith("_tpq.dat")
    return data is not None or split_green is not None or aggregate_green


def series_snapshot(path):
    return {q.name: q.read_bytes() for q in (path / "output").glob("*.dat")
            if is_physical_series(q.name)}


def assert_series_families(snapshot, split_green=False):
    assert snapshot, "physical-series snapshot must be nonempty"
    for family in ("SS", "Norm", "Flct"):
        assert any(re.search(r"(?:^|_){}_".format(family), name)
                   for name in snapshot), (family, sorted(snapshot))
    if split_green:
        assert any(re.search(r"_set[0-9]+step[0-9]+[.]dat$", name)
                   for name in snapshot), sorted(snapshot)


def checkpoint_snapshot(path):
    files = [path / "output/zvo_restart_manifest.dat"]
    files += sorted((path / "output").glob("zvo_restart_set*"))
    assert all(q.is_file() for q in files) and len(files) > 1, files
    return {q.name: q.read_bytes() for q in files}


def append_row_with_boundary(path, name, set_index, step):
    target = path / "output" / name
    lines = target.read_text().splitlines()
    data = next(line.split() for line in reversed(lines)
                if line.strip() and not line.lstrip().startswith("#"))
    data[0], data[1] = str(set_index), str(step)
    target.write_text(target.read_text() + " ".join(data) + "\n")


ref = prepare("ref", 0)
rc.run(ref, "ref")
cont = prepare("cont", 1, interval=4)
rc.run(cont, "interrupted", binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "1:13"})
meta, manifest_rows = rc.manifest(cont)
assert manifest_rows[0]["step"] == 20 and manifest_rows[1]["step"] == 13
rc.set_key(cont, "calcmod.def", "ReStart", 3)
rc.drop_key(cont, "modpara.def", "CheckpointInterval")
rc.run(cont, "continued")
compare_set_files(ref, cont, (0, 1))
assert_ctpq_beta_grid(cont)
if MODE == "ctpq":
    row = rc.manifest(cont)[1][1]
    words, _ = rc.v3_header(rc.restart_files(cont, "zvo", 1, 0,
                                              row["gen"])[0])
    assert word_double(words[36]) > 0
    np.testing.assert_allclose(word_double(words[37]),
                               (row["step"] - 1) / 4.0, rtol=0, atol=1e-15)
    np.testing.assert_allclose(word_double(words[38]), 0.25,
                               rtol=0, atol=1e-15)

agg_ref = prepare("agg_ref", 0, aggregate=True)
rc.run(agg_ref, "ref")
agg = prepare("agg", 1, interval=4, aggregate=True)
rc.run(agg, "interrupted", binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "1:13"})
p = agg / "output/SS_tpq.dat"
p.write_text(p.read_text() + " 1 13 0.1 -1 0 0 0\n 1 14 0.1 -1 0 0 0\n")
rc.set_key(agg, "calcmod.def", "ReStart", 3)
rc.drop_key(agg, "modpara.def", "CheckpointInterval")
rc.run(agg, "continued")
compare_aggregates(agg_ref, agg)
assert_ctpq_beta_grid(agg, aggregate=True)

green_ref = prepare("green_ref", 0, aggregate=True, observations=True)
rc.run(green_ref, "ref")
green = prepare("green", 1, interval=4, aggregate=True, observations=True)
rc.run(green, "interrupted", binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "1:13"})
green_file = green / "output/zvo_cisajs_tpq.dat"
green_file.write_text(green_file.read_text() +
                      " 1 14 0 0 0 0 -1.0 0.0\n")
rc.set_key(green, "calcmod.def", "ReStart", 3)
rc.drop_key(green, "modpara.def", "CheckpointInterval")
rc.run(green, "continued")
observations = sorted_observations(green)
np.testing.assert_allclose(observations, sorted_observations(green_ref),
                           rtol=0, atol=1e-10)
assert observations.size and 0 in observations[:, 1]
if MODE == "ctpq":
    assert 1 in observations[:, 1]
else:
    assert 1 not in observations[:, 1]
assert np.max(observations[:, 1]) >= 2

noop = prepare("noop", 1, aggregate=True)
rc.run(noop, "first")
p = noop / "output/Norm_tpq.dat"
clean = p.read_text()
p.write_text(clean + " 0 20 0.1 1 1\n")
rc.set_key(noop, "calcmod.def", "ReStart", 3)
text = rc.run(noop, "noop")
assert "no remaining steps" in text and p.read_text() == clean

grow = prepare("grow", 1, numave=1, aggregate=True)
rc.run(grow, "one")
before = rc.rows(grow, "SS_tpq.dat")
rc.set_key(grow, "modpara.def", "NumAve", 2)
rc.set_key(grow, "calcmod.def", "ReStart", 2)
text = rc.run(grow, "two")
after = rc.rows(grow, "SS_tpq.dat")
assert "no remaining steps" in text
np.testing.assert_allclose(after[after[:, 0] == 0], before, rtol=0, atol=0)
assert np.sum(after[:, 0] == 1) == before.shape[0]

tiny = prepare("tiny", 1, numave=1, lanczos_max=1)
rc.run(tiny, "first")
assert rc.manifest(tiny)[1][0]["step"] == TINY_STEP
saved = (tiny / "output/SS_rand0.dat").read_bytes()
tiny_checkpoints = checkpoint_snapshot(tiny)
rc.set_key(tiny, "calcmod.def", "ReStart", 3)
text = rc.run(tiny, "noop")
assert "no remaining steps" in text
assert (tiny / "output/SS_rand0.dat").read_bytes() == saved
assert checkpoint_snapshot(tiny) == tiny_checkpoints

# A terminal v3 checkpoint remains a genuine continuation source when the
# absolute iteration limit is raised.
terminal = prepare("terminal_v3", 1, numave=1, lanczos_max=8)
rc.run(terminal, "first")
assert rc.manifest(terminal)[1][0]["step"] == 8
rc.set_key(terminal, "calcmod.def", "ReStart", 3)
rc.set_key(terminal, "modpara.def", "Lanczos_max", 20)
rc.run(terminal, "extended")
ref_terminal = prepare("terminal_ref", 0, numave=1)
rc.run(ref_terminal, "ref")
compare_set_files(ref_terminal, terminal, (0,))

# Lowering the limit makes a completed v3 set a no-op. Recovery removes a
# synthetic tail once; the repeated no-op preserves series and checkpoint bytes.
lowered = prepare("lowered", 1, numave=1)
rc.run(lowered, "first")
checkpoint_paths = ([lowered / "output/zvo_restart_manifest.dat"] +
                    sorted((lowered / "output").glob("zvo_restart_set*")))
checkpoint_bytes = {q.name: q.read_bytes() for q in checkpoint_paths}
for name in ("SS_rand0.dat", "Norm_rand0.dat", "Flct_rand0.dat"):
    q = lowered / "output" / name
    last = q.read_text().splitlines()[-1].split()
    last[-1] = "20"
    q.write_text(q.read_text() + " ".join(last) + "\n")
rc.set_key(lowered, "calcmod.def", "ReStart", 3)
rc.set_key(lowered, "modpara.def", "Lanczos_max", 7)
text = rc.run(lowered, "lowered_noop")
assert "no remaining steps" in text
clean_series = series_snapshot(lowered)
assert {q.name: q.read_bytes() for q in checkpoint_paths} == checkpoint_bytes
rc.run(lowered, "repeated_noop")
assert {q.name: q.read_bytes() for q in checkpoint_paths} == checkpoint_bytes
assert {q.name: q.read_bytes() for q in (lowered / "output").glob("*.dat")
        if q.name in clean_series} == clean_series

leg = prepare("legacy", 1, numave=1, lanczos_max=8)
rc.run(leg, "first")
rc.v3_to_legacy_tmpvec(leg, "zvo", 0)
(leg / "output/zvo_restart_manifest.dat").unlink()
for q in (leg / "output").glob("zvo_restart_set*"):
    q.unlink()
legacy_before = {q.name: q.read_bytes() for q in (leg / "output").glob("tmpvec_set*_rank_*.dat")}
rc.set_key(leg, "calcmod.def", "ReStart", 3)
rc.set_key(leg, "modpara.def", "Lanczos_max", 20)
rc.run(leg, "extend_without_publication")
assert not (leg / "output/zvo_restart_manifest.dat").exists()
assert {q.name: q.read_bytes() for q in (leg / "output").glob("tmpvec_set*_rank_*.dat")} == legacy_before
assert np.max(rc.rows(leg, "SS_rand0.dat")[:, -1]) == 19
rc.run(leg, "resume_from_same_legacy")
ref1 = prepare("ref1", 0, numave=1)
rc.run(ref1, "ref")
for name in ("SS_rand0.dat", "Flct_rand0.dat"):
    np.testing.assert_allclose(rc.rows(leg, name), rc.rows(ref1, name), rtol=0, atol=1e-10)
norm = rc.rows(leg, "Norm_rand0.dat")
norm_ref = rc.rows(ref1, "Norm_rand0.dat")
assert norm.shape == norm_ref.shape
np.testing.assert_allclose(norm[:, [0, 1, 3]], norm_ref[:, [0, 1, 3]], rtol=0, atol=1e-10)
assert np.all(norm[norm[:, -1] >= 8, 2] == 0)


for source_case in (cont, agg):
    for bad_step in (0, 7):
        bad = ROOT / "bad_manifest_{}_{}".format(source_case.name, bad_step)
        shutil.copytree(source_case, bad)
        m = bad / "output/zvo_restart_manifest.dat"
        lines = m.read_text().splitlines()
        lines = [" ".join(line.split()[:7] + [str(bad_step)] + line.split()[8:])
                 if line.startswith("set 1 ") else line for line in lines]
        m.write_text("\n".join(lines) + "\n")
        before = series_snapshot(bad)
        rc.run(bad, "reject_boundary", fail="restart checkpoint")
        assert series_snapshot(bad) == before

for damage, code in (("short", 3), ("step0", 4)):
    bad = ROOT / ("bad_legacy_" + damage)
    shutil.copytree(leg, bad)
    q = sorted((bad / "output").glob("tmpvec_set0_rank_*.dat"))[0]
    data = q.read_bytes()
    q.write_bytes(data[:4] if damage == "short" else struct.pack("=i", 0) + data[4:])
    before = series_snapshot(bad)
    rc.run(bad, "reject_legacy",
           fail="validation failed (code={}".format(code))
    assert series_snapshot(bad) == before

# TPQ FLAGS must be zero. cTPQ additionally rejects negative cumulative beta
# and a delta_tau that differs from 1/LargeValue, all before series recovery.
header_mutations = [("flags", 39, "<Q", 1)]
if MODE == "ctpq":
    header_mutations += [
        ("negative_beta", 37, "<d", -0.25),
        ("delta_tau", 38, "<d", 0.25 * (1.0 + 1e-10)),
    ]
for label, word, fmt, value in header_mutations:
    bad_header = ROOT / ("bad_header_" + label)
    shutil.copytree(cont, bad_header)
    row = rc.manifest(bad_header)[1][1]
    for q in rc.restart_files(bad_header, "zvo", 1, 0, row["gen"]):
        data = bytearray(q.read_bytes())
        struct.pack_into(fmt, data, 8 * word, value)
        q.write_bytes(data)
    before = series_snapshot(bad_header)
    assert_series_families(before)
    rc.run(bad_header, "reject_" + label, fail="restart checkpoint")
    assert series_snapshot(bad_header) == before

np_ = 1
for i, arg in enumerate(rc.MPI[:-1]):
    if arg in ("-np", "-n"):
        np_ = int(rc.MPI[i + 1])

# The same outer-driver discipline is required by both TPQ methods.
solver_name = "cTPQ" if MODE == "ctpq" else "mTPQ"
failure_cases = [
    (MODE + "_output_open", False, solver_name + " output open failed", "output_open"),
    (MODE + "_expec_status", True, solver_name + " driver operation failed", "driver"),
]
if MODE == "tpq":
    failure_cases.append(("tpq_output_close", False,
                          "mTPQ output close failed", "output_close"))
for stage, observations, diagnostic, marker_stage in failure_cases:
    for rank in sorted({0, np_ - 1}):
        bad = prepare("{}_rank{}".format(stage, rank), 0, numave=1,
                      lanczos_max=2, observations=observations)
        text = rc.run(bad, "injected_failure", binary=rc.probe(),
                      fail=diagnostic, timeout=60,
                      env={"HPHI_TEST_FD_FAIL": stage,
                           "HPHI_TEST_FD_FAIL_RANK": str(rank)})
        agreed = {int(value) for value in re.findall(
            r"{} agreed failure stage={} rank=([0-9]+)".format(
                solver_name, marker_stage), text)}
        assert agreed == set(range(np_)), (stage, rank, agreed, np_)
        assert "MPI_ERR" not in text, text

# Actual libc root-open failures and inner observable agreement, including
# publication preservation on a failed ReStart=2 continuation.
import restart_observable_failures as output_failures
output_failures.tpq_failures(ROOT, METHOD)

for rank in sorted({0, np_ - 1}):
    bad = ROOT / "alloc_{}".format(rank)
    shutil.copytree(cont, bad)
    before = series_snapshot(bad)
    rc.run(bad, "alloc", binary=rc.probe(), fail="restart recovery allocation",
           env={"HPHI_TEST_FD_FAIL": "restart_recover_alloc", "HPHI_TEST_FD_FAIL_RANK": str(rank)})
    assert series_snapshot(bad) == before

# Sector-basis continuation at a fixed layout. The same test is launched at
# np1/4/16; both distributed and replicated layouts retain their own saved
# identity throughout each restart chain.
if MODE in ("tpq", "ctpq"):
    import symmetry_general_terms as sector_fixture

    sector_fixture.ROOT = ROOT / "sector_source"
    sector_fixture.ROOT.mkdir()

    def sector_callback(path, model, length, momentum, states, raw, projector):
        del model, length, momentum, states, raw, projector

        def configure(target, restart, aggregate, checkpoint=0,
                      numave=2, lanczos_max=10):
            calc_type = {"tpq": 1, "ctpq": 5}[MODE]
            calc = (target / "calc.def").read_text().replace(
                "CalcType 3", "CalcType {}".format(calc_type))
            (target / "calc.def").write_text(calc +
                "ReStart {}\nInitialVecType -1\nOutputGreenFormat {}\n".format(
                    restart, int(aggregate)))
            mod = (target / "mod.def").read_text().replace(
                "Lanczos_max 400", "Lanczos_max {}".format(lanczos_max)).replace(
                "LargeValue 100", "LargeValue 4")
            (target / "mod.def").write_text(
                mod + "NumAve {}\nExpecInterval {}\nCheckpointInterval {}\n".format(
                    numave, 1 if MODE == "ctpq" else 2, checkpoint))
            (target / "greenone.def").write_text(
                "===============================\nNCisAjs 2\n===============================\n"
                "======== Green functions ======\n===============================\n"
                "0 0 0 0\n1 1 1 1\n")
            sym = target / "sym.def"
            sym.write_text(sym.read_text() + "OneBodyG greenone.def\n")

        def clone(name, restart, aggregate, checkpoint=0,
                  numave=2, lanczos_max=10):
            target = ROOT / name
            shutil.copytree(path, target)
            configure(target, restart, aggregate, checkpoint,
                      numave, lanczos_max)
            return target

        def compare_sector(reference, actual, aggregate, sets=(0, 1)):
            if aggregate:
                for name in ("zvo_SS_tpq.dat", "zvo_Norm_tpq.dat",
                             "zvo_Flct_tpq.dat"):
                    np.testing.assert_allclose(
                        rc.sorted_by_key(rc.rows(actual, name), (0, 1)),
                        rc.sorted_by_key(rc.rows(reference, name), (0, 1)),
                        rtol=0, atol=1e-10)
                np.testing.assert_allclose(
                    sorted_observations(actual), sorted_observations(reference),
                    rtol=0, atol=1e-10)
                assert_ctpq_beta_grid(actual, aggregate=True)
            else:
                for set_index in sets:
                    for family in ("SS", "Norm", "Flct"):
                        name = "zvo_{}_rand{}.dat".format(family, set_index)
                        np.testing.assert_allclose(rc.rows(actual, name),
                                                   rc.rows(reference, name),
                                                   rtol=0, atol=1e-10)
                assert_ctpq_beta_grid(actual, aggregate=False, sets=sets)

        for layout in ("distributed", "replicated"):
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout)

            # Periodic interruption/resume covers both data-output modes and
            # performs real suffix iterations for both samples.
            for aggregate in (False, True):
                stem = "sector_{}_{}".format(layout, "aggregate" if aggregate else "split")
                reference = clone(stem + "_ref", 0, aggregate)
                continued = clone(stem + "_cont", 1, aggregate, checkpoint=3)
                text = rc.run(reference, "ref", env=env, namelist="sym.def")
                assert "Symmetry basis layout: {}".format(layout) in text
                record = dict(line.split("=", 1) for line in
                              (reference / "output/zvo_symmetry_sector.dat").read_text().splitlines())
                assert record["model"] == "Spin" and record["nsite"] == "6"
                assert record["full_dim"] == "20" and record["sector_dim"] == "4"
                assert record["basis_layout"] == layout
                rc.run(continued, "interrupted", binary=rc.probe(), env=dict(
                    env, HPHI_TEST_RESTART_ABORT_AFTER="1:7"), namelist="sym.def")
                rc.set_key(continued, "calc.def", "ReStart", 3)
                rc.set_key(continued, "mod.def", "CheckpointInterval", 0)
                rc.run(continued, "resumed", env=env, namelist="sym.def")
                compare_sector(reference, continued, aggregate)
                if aggregate:
                    if MODE == "ctpq":
                        assert 1 in sorted_observations(continued)[:, 1]
                    else:
                        assert 1 not in sorted_observations(continued)[:, 1]

            # A terminal v3 checkpoint at STEP=6 resumes actual iterations
            # 6..9 when the absolute limit grows to 10.
            terminal_ref = clone("sector_{}_terminal_ref".format(layout), 0, False)
            terminal = clone("sector_{}_terminal".format(layout), 1, False,
                             lanczos_max=6)
            rc.run(terminal_ref, "ref", env=env, namelist="sym.def")
            rc.run(terminal, "first", env=env, namelist="sym.def")
            assert {row["step"] for row in rc.manifest(terminal)[1].values()} == {6}
            rc.set_key(terminal, "calc.def", "ReStart", 3)
            rc.set_key(terminal, "mod.def", "Lanczos_max", 10)
            text = rc.run(terminal, "extended", env=env, namelist="sym.def")
            assert "no remaining steps" not in text
            compare_sector(terminal_ref, terminal, False)
            assert np.max(rc.rows(terminal, "zvo_SS_rand1.dat")[:, -1]) == 9

            # Each method's bootstrap lower bound is a byte-invariant no-op.
            tiny = clone("sector_{}_tiny".format(layout), 1, False,
                         numave=1, lanczos_max=1)
            rc.run(tiny, "first", env=env, namelist="sym.def")
            assert rc.manifest(tiny)[1][0]["step"] == TINY_STEP
            tiny_series = series_snapshot(tiny)
            assert_series_families(tiny_series, split_green=True)
            tiny_checkpoints = checkpoint_snapshot(tiny)
            rc.set_key(tiny, "calc.def", "ReStart", 3)
            text = rc.run(tiny, "noop", env=env, namelist="sym.def")
            assert "no remaining steps" in text
            assert series_snapshot(tiny) == tiny_series
            assert checkpoint_snapshot(tiny) == tiny_checkpoints

            # Lowered-limit recovery removes aggregate physical tails once;
            # the next no-op leaves physical series and checkpoint bytes fixed.
            lowered = clone("sector_{}_lowered".format(layout), 1, True,
                            numave=1)
            rc.run(lowered, "first", env=env, namelist="sym.def")
            lowered_checkpoints = checkpoint_snapshot(lowered)
            for name in ("zvo_SS_tpq.dat", "zvo_Norm_tpq.dat",
                         "zvo_Flct_tpq.dat", "zvo_cisajs_tpq.dat"):
                append_row_with_boundary(lowered, name, 0, 10)
            dirty_series = series_snapshot(lowered)
            assert_series_families(dirty_series)
            rc.set_key(lowered, "calc.def", "ReStart", 3)
            rc.set_key(lowered, "mod.def", "Lanczos_max", 5)
            text = rc.run(lowered, "lowered_noop", env=env, namelist="sym.def")
            assert "no remaining steps" in text
            clean_series = series_snapshot(lowered)
            assert clean_series != dirty_series
            assert checkpoint_snapshot(lowered) == lowered_checkpoints
            rc.run(lowered, "repeated_noop", env=env, namelist="sym.def")
            assert series_snapshot(lowered) == clean_series
            assert checkpoint_snapshot(lowered) == lowered_checkpoints

            # NumAve growth preserves completed set 0 and adds set 1 with the
            # same rows as a two-sample uninterrupted reference.
            growth_ref = clone("sector_{}_growth_ref".format(layout), 0, True)
            growth = clone("sector_{}_growth".format(layout), 1, True,
                           numave=1)
            rc.run(growth_ref, "ref", env=env, namelist="sym.def")
            rc.run(growth, "one", env=env, namelist="sym.def")
            old_row = dict(rc.manifest(growth)[1][0])
            old_set0_files = {q.name: q.read_bytes() for q in
                              rc.restart_files(growth, "zvo", 0, 0,
                                               old_row["gen"])}
            assert old_set0_files
            old_set0 = {name: rc.rows(growth, name)
                        [rc.rows(growth, name)[:, 0] == 0]
                        for name in ("zvo_SS_tpq.dat", "zvo_Norm_tpq.dat",
                                     "zvo_Flct_tpq.dat")}
            assert all(rows.size for rows in old_set0.values())
            rc.set_key(growth, "calc.def", "ReStart", 2)
            rc.set_key(growth, "mod.def", "NumAve", 2)
            text = rc.run(growth, "two", env=env, namelist="sym.def")
            assert "no remaining steps" in text
            rows_after = rc.manifest(growth)[1]
            assert dict(rows_after[0]) == old_row and 1 in rows_after
            assert {q.name: q.read_bytes() for q in
                    rc.restart_files(growth, "zvo", 0, 0,
                                     old_row["gen"])} == old_set0_files
            for name, expected in old_set0.items():
                actual = rc.rows(growth, name)
                np.testing.assert_allclose(actual[actual[:, 0] == 0], expected,
                                           rtol=0, atol=0)
            compare_sector(growth_ref, growth, True)

            # Strict input restart rejects a missing newly-required set before
            # touching prefixed split data or Green observations.
            missing = clone("sector_{}_missing".format(layout), 1, False,
                            numave=1)
            rc.run(missing, "first", env=env, namelist="sym.def")
            missing_series = series_snapshot(missing)
            assert_series_families(missing_series, split_green=True)
            missing_checkpoints = checkpoint_snapshot(missing)
            rc.set_key(missing, "calc.def", "ReStart", 3)
            rc.set_key(missing, "mod.def", "NumAve", 2)
            rc.run(missing, "reject_missing", env=env, namelist="sym.def",
                   fail="requires a checkpoint for every set")
            assert series_snapshot(missing) == missing_series
            assert checkpoint_snapshot(missing) == missing_checkpoints

            # Independent set publication may leave a valid mixed-layout
            # manifest. Exercise BOTH directions with two existing sets.
            other = "replicated" if layout == "distributed" else "distributed"
            target_env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=other)
            switch = clone("switch_" + layout, 1, True, lanczos_max=6)
            rc.run(switch, "initial", env=env, namelist="sym.def")
            rc.set_key(switch, "calc.def", "ReStart", 2)
            rc.set_key(switch, "mod.def", "Lanczos_max", 10)
            rc.set_key(switch, "mod.def", "CheckpointInterval", 3)
            text = rc.run(switch, "switch_set0", binary=rc.probe(),
                          env=dict(target_env, HPHI_TEST_RESTART_ABORT_AFTER="0:7"),
                          namelist="sym.def")
            rows = rc.manifest(switch)[1]
            assert rows[0]["step"] == 7 and rows[1]["step"] == 6, rows
            for set_index, expected_layout in ((0, other), (1, layout)):
                headers = [rc.v3_header(q)[0] for q in rc.restart_files(
                    switch, "zvo", set_index, 0, rows[set_index]["gen"])]
                assert len(headers) == np_
                assert {h[13] for h in headers} == {
                    1 if expected_layout == "distributed" else 0}
                if np_ > 4:
                    assert any(h[15] == 0 for h in headers), headers
            for restart in (2, 3):
                resumed = ROOT / ("switch_{}_mode{}".format(layout, restart))
                shutil.copytree(switch, resumed)
                before = checkpoint_snapshot(resumed)
                rc.set_key(resumed, "calc.def", "ReStart", restart)
                if restart == 3:
                    rc.set_key(resumed, "mod.def", "CheckpointInterval", 0)
                rc.run(resumed, "continue", env=target_env, namelist="sym.def")
                compare_sector(growth_ref, resumed, True)
                if restart == 3:
                    assert checkpoint_snapshot(resumed) == before
                else:
                    assert {r["step"] for r in rc.manifest(resumed)[1].values()} == {10}

        # Logical identity must match the runtime basis even though metadata
        # layout is compatible. Invalid layout and within-set disagreement
        # remain rejected before any output recovery.
        for damage in ("identity", "invalid_layout", "within_set_layout"):
            if damage == "within_set_layout" and np_ == 1:
                continue
            bad = clone("sector_bad_" + damage, 1, False)
            env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT="replicated")
            rc.run(bad, "first", env=env, namelist="sym.def")
            row = rc.manifest(bad)[1][1]
            files = rc.restart_files(bad, "zvo", 1, 0, row["gen"])
            if damage == "within_set_layout":
                files = files[-1:]
            for q in files:
                data = bytearray(q.read_bytes())
                word = 16 if damage == "identity" else 13
                value = struct.unpack_from("<Q", data, 8 * word)[0] ^ 1
                if damage == "invalid_layout":
                    value = 2
                struct.pack_into("<Q", data, 8 * word, value)
                q.write_bytes(data)
            before = series_snapshot(bad)
            assert_series_families(before, split_green=True)
            rc.set_key(bad, "calc.def", "ReStart", 3)
            rc.run(bad, "reject_" + damage, env=env, namelist="sym.def",
                   fail="restart checkpoint")
            assert series_snapshot(bad) == before

    sector_fixture.prepare("Spin", 6, sector_test=sector_callback,
                           sector_momenta=(0,))

print("ok")

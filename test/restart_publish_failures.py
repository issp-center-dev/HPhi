"""Checkpoint publication failures preserve the last complete unit."""
import re
import shutil

import restart_common as rc

ROOT = rc.fresh_dir("restart_publish_failures")
STAN = ('L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n'
        '2Sz=0\nexct=2\nLanczos_max=40\noutputmode="None"\n')
NP = rc.mpi_size()


def prepare(name):
    path = ROOT / name
    rc.generate(path, STAN)
    rc.set_key(path, "calcmod.def", "ReStart", 2)
    rc.set_key(path, "modpara.def", "CheckpointInterval", 5)
    return path


def snapshot(path):
    return {p.name: p.read_bytes()
            for p in (path / "output").glob("zvo_restart_*")}


def assert_published_files(path):
    found = rc.manifest(path)
    assert found is not None, path
    _, rows = found
    expected = {"zvo_restart_manifest.dat"}
    for set_index, row in rows.items():
        for comp in range(row["ncomp"]):
            expected.update(
                "zvo_restart_set{}_comp{}_gen{}_rank_{}.dat".format(
                    set_index, comp, row["gen"], rank)
                for rank in range(row["ranks"]))
    actual = {p.name for p in (path / "output").glob("zvo_restart_*")}
    assert actual == expected, (sorted(actual), sorted(expected))


cases = [
    ("restart_write", "0", "0"),
    ("restart_rename", "0", "0"),
    ("restart_manifest", "0", "0"),
    ("restart_write", "0", "3"),       # state 1, direction component
    ("restart_rename", "0", "2"),      # state 1, state component
]
if NP > 1:
    cases += [
        ("restart_write", str(NP - 1), "1"),
        ("restart_rename", str(NP - 1), "3"),
        ("restart_manifest", str(NP - 1), "0"),
    ]

agreement_stage = {
    "restart_write": "write",
    "restart_rename": "rename",
    "restart_manifest": "manifest update",
}
for stage, rank, skip in cases:
    path = prepare("{}_{}_{}".format(stage, rank, skip))
    rc.run(path, "first", binary=rc.probe(),
           env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:6"})
    before = snapshot(path)
    text = rc.run(
        path, stage, binary=rc.probe(),
        env={"HPHI_TEST_FD_FAIL": stage,
             "HPHI_TEST_FD_FAIL_RANK": rank,
             "HPHI_TEST_FD_FAIL_SKIP": skip},
        fail=("restart manifest" if stage == "restart_manifest"
              else "restart checkpoint"), timeout=120)
    agreed = {int(value) for value in re.findall(
        r"restart publication agreed failure stage={} rank=([0-9]+)".format(
            re.escape(agreement_stage[stage])), text)}
    assert agreed == set(range(NP)), (stage, rank, skip, agreed, NP)
    assert snapshot(path) == before, (stage, rank, skip)
    rc.run(path, "resume")
    assert_published_files(path)

# Seed the generation not referenced by the manifest and a stale .part.  The
# manifest generation remains the source of truth; the next checkpoint safely
# replaces the stale generation and removes the superseded published files.
stale = prepare("stale_unreferenced")
rc.run(stale, "first", binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:6"})
old_rows = rc.manifest(stale)[1]
old_gen = old_rows[0]["gen"]
assert {row["gen"] for row in old_rows.values()} == {old_gen}
new_gen = 1 - old_gen
for old in list((stale / "output").glob(
        "zvo_restart_set*_comp*_gen{}_rank_*.dat".format(old_gen))):
    replacement = old.with_name(old.name.replace(
        "_gen{}_".format(old_gen), "_gen{}_".format(new_gen)))
    replacement.write_bytes(b"stale unreferenced generation\n")
part = stale / "output" / ("zvo_restart_set0_comp0_gen{}_rank_0.dat.part".format(new_gen))
part.write_bytes(b"stale partial publication\n")
rc.run(stale, "continued", binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:11"})
stale_rows = rc.manifest(stale)[1]
assert {row["gen"] for row in stale_rows.values()} == {new_gen}
assert {row["step"] for row in stale_rows.values()} == {11}
assert_published_files(stale)

# A normal later publication retains only the manifest-referenced generation.
leftover = prepare("leftover")
rc.run(leftover, "first", binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:11"})
assert_published_files(leftover)

# Fresh ReStart=1 in a directory holding another method's manifest replaces it.
tpq = ROOT / "fresh_tpq"
rc.generate(tpq, STAN.replace('method="CG"', 'method="TPQ"').replace(
    "exct=2\n", "NumAve=1\n"))
(tpq / "output").mkdir(exist_ok=True)
for p in (leftover / "output").glob("zvo_restart_*"):
    shutil.copyfile(p, tpq / "output" / p.name)
rc.set_key(tpq, "calcmod.def", "ReStart", 1)
rc.run(tpq, "fresh")
meta, rows = rc.manifest(tpq)
assert meta["method"] == "TPQ" and set(rows) == {0}, (meta, rows)
assert_published_files(tpq)

# Fresh same-method TPQ with reduced NumAve discards old unprocessed rows at
# the first publication rather than waiting for the terminal set.
old_same_tpq = ROOT / "old_same_tpq"
rc.generate(old_same_tpq, STAN.replace('method="CG"', 'method="TPQ"').replace(
    "exct=2\n", "NumAve=2\n"))
rc.set_key(old_same_tpq, "calcmod.def", "ReStart", 1)
rc.run(old_same_tpq, "save")
assert set(rc.manifest(old_same_tpq)[1]) == {0, 1}
reduced = ROOT / "fresh_same_tpq_reduced"
rc.generate(reduced, STAN.replace('method="CG"', 'method="TPQ"').replace(
    "exct=2\n", "NumAve=1\n"))
(reduced / "output").mkdir(exist_ok=True)
for p in (old_same_tpq / "output").glob("zvo_restart_*"):
    shutil.copyfile(p, reduced / "output" / p.name)
rc.set_key(reduced, "calcmod.def", "ReStart", 1)
rc.set_key(reduced, "modpara.def", "CheckpointInterval", 5)
text = rc.run(reduced, "first_publication", binary=rc.probe(),
              env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:6"})
rows = rc.manifest(reduced)[1]
assert set(rows) == {0} and rows[0]["step"] == 6, rows
assert "Test hook: stopping after checkpoint publication." in text
assert not list((reduced / "output").glob("zvo_restart_set1_*"))
assert_published_files(reduced)

# Safe fresh TPQ -> CG: every injected stage preserves every old payload byte.
old_tpq = ROOT / "old_tpq"
rc.generate(old_tpq, STAN.replace('method="CG"', 'method="TPQ"').replace(
    "exct=2\n", "NumAve=2\n"))
rc.set_key(old_tpq, "calcmod.def", "ReStart", 1)
rc.run(old_tpq, "save")
assert {r["gen"] for r in rc.manifest(old_tpq)[1].values()} == {0}
for stage in ("restart_write", "restart_rename", "restart_manifest"):
    case = prepare("fresh_" + stage)
    rc.set_key(case, "calcmod.def", "ReStart", 1)
    (case / "output").mkdir(exist_ok=True)
    for q in (old_tpq / "output").glob("zvo_restart_*"):
        shutil.copyfile(q, case / "output" / q.name)
    before = snapshot(case)
    text = rc.run(
        case, "fail", binary=rc.probe(),
        env={"HPHI_TEST_FD_FAIL": stage, "HPHI_TEST_FD_FAIL_RANK": "0"},
        fail=("restart manifest" if stage == "restart_manifest"
              else "restart checkpoint"), timeout=120)
    agreed = {int(value) for value in re.findall(
        r"restart publication agreed failure stage={} rank=([0-9]+)".format(
            re.escape(agreement_stage[stage])), text)}
    assert agreed == set(range(NP)), (stage, agreed, NP)
    assert snapshot(case) == before

# A mixed-parity old TPQ unit occupies both possible generations, so a fresh
# two-state CG unit must refuse publication without changing any old byte.
mixed = ROOT / "mixed_old_tpq"
shutil.copytree(old_tpq, mixed)
for q in list((mixed / "output").glob("zvo_restart_set1_*_gen0_*")):
    q.rename(q.with_name(q.name.replace("_gen0_", "_gen1_")))
manifest = mixed / "output/zvo_restart_manifest.dat"
manifest.write_text("\n".join(
    line.replace(" gen 0 ", " gen 1 ") if line.startswith("set 1 ") else line
    for line in manifest.read_text().splitlines()) + "\n")
rc.set_key(mixed, "calcmod.def", "ReStart", 3)
rc.run(mixed, "validate_old")
case = prepare("fresh_mixed_refused")
rc.set_key(case, "calcmod.def", "ReStart", 1)
(case / "output").mkdir(exist_ok=True)
for q in (mixed / "output").glob("zvo_restart_*"):
    shutil.copyfile(q, case / "output" / q.name)
before = snapshot(case)
rc.run(case, "refused", fail="no collision-free generation")
assert snapshot(case) == before

print("restart publication failures: ok (np={})".format(NP), flush=True)

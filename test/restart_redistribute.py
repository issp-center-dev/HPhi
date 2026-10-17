"""v3 checkpoints continue with a different MPI rank count and sector layout."""
import os
import shutil
from pathlib import Path

import numpy as np

import restart_common as rc
import symmetry_general_terms as fixture


ROOT = rc.fresh_dir("restart_redistribute")
CG = ('L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\n'
      'exct=2\nLanczos_max=40\nLanczosEps=10\noutputmode="None"\n')
TPQ = ('L=8\nmodel="Spin"\nmethod="TPQ"\nlattice="chain"\nJ=1\n2Sz=0\n'
       'NumAve=2\nLanczos_max=20\nLargeValue=4\nExpecInterval=1000\n'
       'outputmode="None"\n')
ENV = {"OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "1")}


def run_np(path, label, ranks, **kwargs):
    saved = rc.MPI
    rc.MPI = rc.launcher(ranks)
    env = dict(ENV)
    env.update(kwargs.pop("env", {}))
    try:
        return rc.run(path, label, env=env, **kwargs)
    finally:
        rc.MPI = saved


def steps(path):
    return np.loadtxt(path / "output/zvo_Lanczos_Step.dat", skiprows=1, ndmin=2)


def same_history(actual, expected):
    assert actual.shape == expected.shape, (actual.shape, expected.shape)
    np.testing.assert_allclose(actual[:, 3:], expected[:, 3:], rtol=0, atol=1e-9)
    np.testing.assert_allclose(actual[:, 1], expected[:, 1], rtol=1e-6, atol=1e-12)


# The random initial vector depends on rank/thread count.  Compare every
# continuation with an uninterrupted run made with the writer's rank count.
references = {}
for source_ranks, target_ranks in ((4, 1), (4, 16), (16, 4), (1, 4)):
    if source_ranks not in references:
        ref = ROOT / "cg_ref_{}".format(source_ranks)
        rc.generate(ref, CG)
        run_np(ref, "ref", source_ranks)
        references[source_ranks] = steps(ref)
    case = ROOT / "cg_{}_to_{}".format(source_ranks, target_ranks)
    rc.generate(case, CG)
    rc.set_key(case, "calcmod.def", "ReStart", 1)
    rc.set_key(case, "modpara.def", "CheckpointInterval", 5)
    run_np(case, "save", source_ranks, binary=rc.probe(),
           env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:11"})
    row = rc.manifest(case)[1][0]
    assert row["ranks"] == source_ranks, row
    manifest_bytes = (case / "output/zvo_restart_manifest.dat").read_bytes()
    rc.set_key(case, "calcmod.def", "ReStart", 3)
    rc.drop_key(case, "modpara.def", "CheckpointInterval")
    text = run_np(case, "resume", target_ranks)
    same_history(steps(case), references[source_ranks])
    assert (case / "output/zvo_restart_manifest.dat").read_bytes() == manifest_bytes
    assert "payload files opened" in text


# A partial TPQ republication may leave rows written by different rank counts.
tpq_ref = ROOT / "tpq_ref_4"
rc.generate(tpq_ref, TPQ)
run_np(tpq_ref, "ref", 4)
tpq = ROOT / "tpq_4_to_2"
rc.generate(tpq, TPQ)
rc.set_key(tpq, "calcmod.def", "ReStart", 1)
rc.set_key(tpq, "modpara.def", "CheckpointInterval", 4)
run_np(tpq, "save", 4, binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "1:13"})
rc.set_key(tpq, "calcmod.def", "ReStart", 2)
run_np(tpq, "partial", 2, binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "1:17"})
rows = rc.manifest(tpq)[1]
assert rows[0]["ranks"] == 4 and rows[1]["ranks"] == 2, rows
run_np(tpq, "finish", 2)
for sample in range(2):
    np.testing.assert_allclose(rc.rows(tpq, "SS_rand{}.dat".format(sample)),
                               rc.rows(tpq_ref, "SS_rand{}.dat".format(sample)),
                               rtol=0, atol=1e-10)


# Sector checkpoint: change both rank count and basis layout, then read with
# more ranks than vector elements so that empty owners participate.
def sector_scenario(path, model, length, momentum, states, raw, projector):
    rc.set_key(path, "calc.def", "ReStart", 1)
    rc.set_key(path, "mod.def", "CheckpointInterval", 1)
    rc.set_key(path, "mod.def", "exct", 1)
    ref = Path(str(path) + "_ref")
    shutil.copytree(path, ref)
    rc.set_key(ref, "calc.def", "ReStart", 0)
    rc.drop_key(ref, "mod.def", "CheckpointInterval")
    run_np(ref, "ref", 3, namelist="sym.def",
           env={"HPHI_SYMMETRY_BASIS_LAYOUT": "distributed"})
    saved_text = run_np(path, "save", 3, namelist="sym.def", binary=rc.probe(),
                        env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:2",
                             "HPHI_SYMMETRY_BASIS_LAYOUT": "distributed"})
    assert "Test hook: stopping after checkpoint publication" in saved_text
    row = rc.manifest(path)[1][0]
    assert row["step"] == 2, row
    words, _ = rc.v3_header(rc.restart_files(path, "zvo", 0, 0, row["gen"])[0])
    assert words[39] == 0
    assert steps(ref)[-1, 0] >= 2 and steps(path)[-1, 0] == 1
    rc.set_key(path, "calc.def", "ReStart", 3)
    rc.drop_key(path, "mod.def", "CheckpointInterval")
    resumed = run_np(path, "resume6", 6, namelist="sym.def",
                     env={"HPHI_SYMMETRY_BASIS_LAYOUT": "replicated"})
    assert "no remaining steps" not in resumed and steps(path)[-1, 0] > 1
    same_history(steps(path), steps(ref))
    run_np(path, "resume16", 16, namelist="sym.def",
           env={"HPHI_SYMMETRY_BASIS_LAYOUT": "distributed"})
    same_history(steps(path), steps(ref))


fixture.ROOT = ROOT / "sector"
fixture.ROOT.mkdir()
fixture.prepare("Spin", 6, sector_test=sector_scenario, sector_momenta=(0,))
print("v3 restart rank redistribution PASS")

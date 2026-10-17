"""Legacy raw rank redistribution and remaining v3 method coverage."""
import os
from pathlib import Path
import shutil
import struct

import numpy as np

import restart_common as rc
import restart_te_contract as te_fixture


ROOT = rc.fresh_dir("restart_legacy_sources")
CG = ('L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\n'
      'exct=1\nLanczos_max=100\noutputmode="None"\n')
TPQ = ('L=8\nmodel="Spin"\nmethod="TPQ"\nlattice="chain"\nJ=1\n2Sz=0\n'
       'NumAve=1\nLanczos_max=20\nLargeValue=4\nExpecInterval=1000\n'
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


def energy(path):
    return [float(line.split()[1])
            for line in (path / "output/zvo_energy.dat").read_text().splitlines()
            if line.lstrip().startswith("Energy")]


def build_tjgc(path):
    """Two-site expert tJGC fixture with one empty rank at np4."""
    path.mkdir()
    (path / "namelist.def").write_text(
        "ModPara modpara.def\nLocSpin locspn.def\nTrans trans.def\nCalcMod calcmod.def\n")
    (path / "calcmod.def").write_text(
        "CalcType 3\nCalcModel 10\nReStart 0\nInputEigenVec 0\nOutputEigenVec 0\n")
    (path / "modpara.def").write_text(
        "--------------------\nModel_Parameters 0\n--------------------\n"
        "HPhi_Cal_Parameters\n--------------------\nCDataFileHead zvo\n"
        "CParaFileHead zqp\n--------------------\nNsite 2\nLanczos_max 100\n"
        "initial_iv 1\nexct 1\nLanczosEps 12\nLanczosTarget 2\nLargeValue 12\n")
    (path / "locspn.def").write_text(
        "====\nNlocalSpin 0\n====\n====\n====\n0 0\n1 0\n")
    (path / "trans.def").write_text(
        "====\nNTransfer 4\n====\n====\n====\n"
        "0 0 1 0 1 0\n1 0 0 0 1 0\n0 1 1 1 1 0\n1 1 0 1 1 0\n")


# Raw OutputEigenVec written by four ranks is imported by smaller and larger
# runs only when its source rank count is explicit.
case = ROOT / "eigenvec"
rc.generate(case, CG)
rc.set_key(case, "calcmod.def", "OutputEigenVec", 1)
run_np(case, "write", 4)
energy4 = energy(case)
rc.set_key(case, "calcmod.def", "OutputEigenVec", 0)
rc.set_key(case, "calcmod.def", "InputEigenVec", 1)
run_np(case, "missing_key", 1, fail="RestartSourceRanks is required")
rc.set_key(case, "modpara.def", "RestartSourceRanks", 4)
for target_ranks in (1, 16):
    text = run_np(case, "import{}".format(target_ranks), target_ranks)
    assert "Legacy raw files carry no identity words" in text
    np.testing.assert_allclose(energy(case), energy4, rtol=0, atol=1e-10)
missing = case / "output/zvo_eigenvec_0_rank_3.dat"
saved_missing = missing.read_bytes()
missing.unlink()
run_np(case, "missing_rank", 1, fail="rank_3")
missing.write_bytes(saved_missing)


# An explicitly declared empty trailing writer rank is still mandatory.  This
# tJGC L=2/np4 partition has local dimensions [3, 3, 3, 0].
empty = ROOT / "empty_rank"
build_tjgc(empty)
rc.set_key(empty, "calcmod.def", "OutputEigenVec", 1)
run_np(empty, "write", 4)
empty_files = sorted((empty / "output").glob("zvo_eigenvec_0_rank_*.dat"))
assert [struct.unpack_from("=Q", p.read_bytes(), 4)[0] for p in empty_files] == [3, 3, 3, 0]
empty_files[-1].unlink()
rc.set_key(empty, "calcmod.def", "OutputEigenVec", 0)
rc.set_key(empty, "calcmod.def", "InputEigenVec", 1)
rc.set_key(empty, "modpara.def", "RestartSourceRanks", 4)
run_np(empty, "missing_empty_rank", 1, fail="rank_3")


# Legacy tmpvec continuation (mTPQ) np4 -> np2, compared with an uninterrupted
# four-rank reference whose random seed and prefix are identical.
reference = ROOT / "tpq_ref"
rc.generate(reference, TPQ)
run_np(reference, "ref", 4)
tpq = ROOT / "tpq"
rc.generate(tpq, TPQ)
rc.set_key(tpq, "calcmod.def", "ReStart", 1)
rc.set_key(tpq, "modpara.def", "Lanczos_max", 8)
run_np(tpq, "first", 4)
rc.v3_to_legacy_tmpvec(tpq, "zvo", 0)
(tpq / "output/zvo_restart_manifest.dat").unlink()
for checkpoint in (tpq / "output").glob("zvo_restart_set*"):
    checkpoint.unlink()
rc.set_key(tpq, "calcmod.def", "ReStart", 3)
rc.set_key(tpq, "modpara.def", "Lanczos_max", 20)
rc.set_key(tpq, "modpara.def", "RestartSourceRanks", 4)
run_np(tpq, "resume", 2)
expected = rc.rows(reference, "SS_rand0.dat")
actual = rc.rows(tpq, "SS_rand0.dat")
np.testing.assert_allclose(actual[actual[:, -1] >= 8][:, :2],
                           expected[expected[:, -1] >= 8][:, :2],
                           rtol=0, atol=1e-10)


# cTPQ v3 continuation completes the TPQ-family rank-redistribution matrix.
ctpq_input = TPQ.replace('method="TPQ"', 'method="cTPQ"')
ctpq_reference = ROOT / "ctpq_ref"
rc.generate(ctpq_reference, ctpq_input)
run_np(ctpq_reference, "ref", 4)
ctpq = ROOT / "ctpq"
rc.generate(ctpq, ctpq_input)
rc.set_key(ctpq, "calcmod.def", "ReStart", 1)
rc.set_key(ctpq, "modpara.def", "CheckpointInterval", 4)
run_np(ctpq, "save", 4, binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:9"})
rc.set_key(ctpq, "calcmod.def", "ReStart", 3)
rc.drop_key(ctpq, "modpara.def", "CheckpointInterval")
run_np(ctpq, "resume", 2)
np.testing.assert_allclose(rc.rows(ctpq, "SS_rand0.dat"),
                           rc.rows(ctpq_reference, "SS_rand0.dat"),
                           rtol=0, atol=1e-10)


# The seed CG and the TE reference/save are deliberately fixed at np4 even
# when this outer CTest instance uses np16.
times = [0, .03, .09, .09, .13, .2, .25, .3]
te_reference = te_fixture.build_time_table_case(ROOT / "te_ref", times, writer_ranks=4)
run_np(te_reference, "ref", 4)
te = te_fixture.build_time_table_case(ROOT / "te", times, writer_ranks=4)
rc.set_key(te, "calcmod.def", "ReStart", 1)
run_np(te, "save", 4, binary=rc.probe(),
       env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:6"})
rc.set_key(te, "calcmod.def", "ReStart", 3)
rc.set_key(te, "calcmod.def", "InputEigenVec", 0)
run_np(te, "resume", 1)
for name in ("SS.dat", "Norm.dat", "Flct.dat"):
    np.testing.assert_allclose(rc.rows(te, name), rc.rows(te_reference, name),
                               rtol=0, atol=5e-12, err_msg=name)

print("legacy raw and cTPQ/TE rank redistribution PASS")

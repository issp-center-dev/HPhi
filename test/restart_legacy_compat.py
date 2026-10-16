"""Legacy restart and unchanged raw/spectrum format compatibility."""
import shutil
import struct

import numpy as np

import restart_common as rc


ROOT = rc.fresh_dir("restart_legacy_compat")


def clone(source, name):
    target = ROOT / name
    shutil.copytree(source, target)
    return target


def remove_v3(path):
    manifest = path / "output/zvo_restart_manifest.dat"
    if manifest.exists():
        manifest.unlink()
    for checkpoint in (path / "output").glob("zvo_restart_set*"):
        checkpoint.unlink()


def assert_raw12(path, pattern):
    files = sorted((path / "output").glob(pattern),
                   key=lambda p: int(p.stem.split("rank_")[1]))
    assert files, (path, pattern)
    for checkpoint in files:
        data = checkpoint.read_bytes()
        step, dim = struct.unpack_from("=iQ", data)
        assert step > 0 and len(data) == 12 + 16 * (dim + 1), (
            checkpoint, step, dim, len(data))
        assert data[:8] != struct.pack("<Q", rc.MAGIC), checkpoint


def cg_legacy_additional_iterations():
    case = ROOT / "cg_legacy"
    rc.generate(case,
        'L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n'
        '2Sz=0\nexct=1\nLanczos_max=4\nLanczosEps=10\n'
        'EigenVecIO="out"\noutputmode="None"\n')
    rc.set_key(case, "calcmod.def", "ReStart", 1)
    rc.run(case, "saved")
    row = rc.manifest(case)[1][0]
    assert row["step"] == 5, row
    rc.v3_to_legacy_tmpvec(case, "zvo", 0)
    legacy = sorted((case / "output").glob("tmpvec_set0_rank_*.dat"))
    assert legacy and {struct.unpack_from("=i", q.read_bytes())[0]
                       for q in legacy} == {5}
    step_file = case / "output/zvo_Lanczos_Step.dat"
    step_file.write_text(step_file.read_text() +
                         "      999  1.0  1.0  -999.0\n")
    remove_v3(case)
    rc.set_key(case, "calcmod.def", "ReStart", 3)
    rc.set_key(case, "calcmod.def", "OutputEigenVec", 1)
    rc.set_key(case, "modpara.def", "Lanczos_max", 40)
    rc.run(case, "additional")
    history = np.loadtxt(case / "output/zvo_Lanczos_Step.dat",
                         skiprows=1, ndmin=2)
    assert history.shape[0] > 5
    np.testing.assert_array_equal(history[:, 0],
                                  np.arange(1, int(history[-1, 0]) + 1))
    assert 999 not in history[:, 0]
    raw = sorted((case / "output").glob("zvo_eigenvec_0_rank_*.dat"))
    assert raw and {struct.unpack_from("=i", q.read_bytes())[0]
                    for q in raw} == {int(history[-1, 0])}
    assert_raw12(case, "zvo_eigenvec_0_rank_*.dat")


def prepare_tpq(name, method, lanczos_max, restart, numave=2):
    path = ROOT / name
    rc.generate(path,
        'L=8\nmodel="Spin"\nmethod="{}"\nlattice="chain"\nJ=1\n'
        '2Sz=0\nNumAve={}\nLanczos_max={}\nLargeValue=4\n'
        'ExpecInterval=1000\nRestart="{}"\noutputmode="None"\n'.format(
            method, numave, lanczos_max, restart))
    return path


def append_tpq_tail(path, set_index, steps):
    for family in ("SS", "Norm", "Flct"):
        target = path / "output/{}_rand{}.dat".format(family, set_index)
        lines = target.read_text().splitlines()
        row = next(line.split() for line in reversed(lines)
                   if line.strip() and not line.lstrip().startswith("#"))
        extra = []
        for step in steps:
            marked = list(row)
            marked[0] = "9876.5"
            marked[-1] = str(step)
            extra.append(" ".join(marked))
        target.write_text(target.read_text() + "\n".join(extra) + "\n")


def compare_tpq(reference, actual, set_index, legacy_suffix=False):
    for family in ("SS", "Flct"):
        np.testing.assert_allclose(
            rc.rows(actual, "{}_rand{}.dat".format(family, set_index)),
            rc.rows(reference, "{}_rand{}.dat".format(family, set_index)),
            rtol=0, atol=1e-10, err_msg="{} set {}".format(family, set_index))
    expected = rc.rows(reference, "Norm_rand{}.dat".format(set_index))
    found = rc.rows(actual, "Norm_rand{}.dat".format(set_index))
    if legacy_suffix:
        np.testing.assert_allclose(found[:, [0, 1, 3]],
                                   expected[:, [0, 1, 3]],
                                   rtol=0, atol=1e-10)
        assert np.all(found[found[:, -1] >= 6, 2] == 0.0), found
    else:
        np.testing.assert_allclose(found, expected, rtol=0, atol=1e-10)


def tpq_legacy_compat(method):
    tag = method.lower()
    reference = prepare_tpq(tag + "_reference", method, 10, "none")
    rc.run(reference, "reference")
    source = prepare_tpq(tag + "_source", method, 6, "restart_out")
    rc.run(source, "source")
    for set_index in (0, 1):
        rc.v3_to_legacy_tmpvec(source, "zvo", set_index)

    mixed = clone(source, tag + "_mixed_v3_legacy")
    manifest = mixed / "output/zvo_restart_manifest.dat"
    manifest.write_text("".join(line for line in manifest.read_text().splitlines(True)
                                if not line.startswith("set 1 ")))
    append_tpq_tail(mixed, 0, (6, 7))
    append_tpq_tail(mixed, 1, (6, 7))
    rc.set_key(mixed, "calcmod.def", "ReStart", 3)
    rc.set_key(mixed, "modpara.def", "Lanczos_max", 10)
    rc.run(mixed, "mixed")
    assert set(rc.manifest(mixed)[1]) == {0}
    compare_tpq(reference, mixed, 0)
    compare_tpq(reference, mixed, 1, legacy_suffix=True)
    if method == "cTPQ":
        for set_index in (0, 1):
            table = rc.rows(mixed, "SS_rand{}.dat".format(set_index))
            np.testing.assert_allclose(table[:, 0], table[:, -1] / 4.0,
                                       rtol=0, atol=1e-12)

    # Simulate work after the last legacy STEP that stopped before any v3
    # publication.  The next attempt must trim that tail and regenerate it.
    tails = clone(source, tag + "_legacy_tails")
    remove_v3(tails)
    for q in (tails / "output").glob("tmpvec_set1_rank_*.dat"):
        q.unlink()
    rc.set_key(tails, "modpara.def", "NumAve", 1)
    rc.set_key(tails, "calcmod.def", "ReStart", 3)
    rc.set_key(tails, "modpara.def", "Lanczos_max", 8)
    rc.run(tails, "unpublished_extension")
    assert not (tails / "output/zvo_restart_manifest.dat").exists()
    append_tpq_tail(tails, 0, (7,))
    assert np.any(rc.rows(tails, "SS_rand0.dat")[:, 0] == 9876.5)
    rc.set_key(tails, "modpara.def", "Lanczos_max", 10)
    rc.run(tails, "resume_same_legacy")
    assert not np.any(rc.rows(tails, "SS_rand0.dat")[:, 0] == 9876.5)
    compare_tpq(reference, tails, 0, legacy_suffix=True)


def corrupt_v3_is_authoritative():
    case = prepare_tpq("corrupt_v3_with_legacy", "TPQ", 6,
                       "restart_out", numave=1)
    rc.run(case, "source")
    rc.v3_to_legacy_tmpvec(case, "zvo", 0)
    row = rc.manifest(case)[1][0]
    target = rc.restart_files(case, "zvo", 0, 0, row["gen"])[0]
    data = bytearray(target.read_bytes())
    struct.pack_into("<Q", data, 8, 99)
    target.write_bytes(data)
    legacy_before = {q.name: q.read_bytes() for q in
                     (case / "output").glob("tmpvec_set0_rank_*.dat")}
    rc.set_key(case, "calcmod.def", "ReStart", 3)
    rc.set_key(case, "modpara.def", "Lanczos_max", 10)
    text = rc.run(case, "reject", fail="restart checkpoint header scan failed")
    assert "legacy recovery payload" not in text
    assert {q.name: q.read_bytes() for q in
            (case / "output").glob("tmpvec_set0_rank_*.dat")} == legacy_before


def spectrum_input(calc_spec, lanczos_max):
    return (
        'L=8\nmodel="Spin"\nmethod="Lanczos"\nlattice="chain"\nJ=1\n'
        '2Sz=0\nexct=1\nLanczos_max={}\nLanczosEps=16\n'
        'EigenVecIO="out"\noutputmode="None"\nSpectrumQL=0.5\n'
        'NOmega=5\nOmegaIm=1.0\nCalcSpec="{}"\n'
        'SpectrumType="SzSz"\n').format(lanczos_max, calc_spec)


def seed_spectrum(path, ground):
    (path / "output").mkdir(exist_ok=True)
    for q in (ground / "output").glob("zvo_eigenvec_0_rank_*.dat"):
        shutil.copyfile(q, path / "output" / q.name)
    shutil.copyfile(ground / "output/zvo_energy.dat",
                    path / "output/zvo_energy.dat")


def spectrum_legacy_restart():
    ground = ROOT / "spectrum_ground"
    rc.generate(ground,
        'L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n'
        '2Sz=0\nexct=1\nLanczos_max=80\nLanczosEps=12\n'
        'EigenVecIO="out"\noutputmode="None"\nSpectrumQL=0.5\n'
        'NOmega=5\nOmegaIm=1.0\n')
    rc.run(ground, "ground")
    assert_raw12(ground, "zvo_eigenvec_0_rank_*.dat")

    normal = ROOT / "spectrum_normal"
    rc.generate(normal, spectrum_input("Normal", 20))
    seed_spectrum(normal, ground)
    rc.run(normal, "normal")

    saved = ROOT / "spectrum_saved2"
    rc.generate(saved, spectrum_input("restart_out", 2))
    seed_spectrum(saved, ground)
    rc.run(saved, "saved")
    expected = rc.rows(normal, "zvo_DynamicalGreen.dat")
    coarse = rc.rows(saved, "zvo_DynamicalGreen.dat")
    assert np.max(np.abs(coarse - expected)) > 1e-7

    for mode, value in (("restart_in", 4), ("restart", 5)):
        continued = clone(saved, "spectrum_calcspec{}".format(value))
        rc.set_key(continued, "calcmod.def", "CalcSpec", value)
        rc.set_key(continued, "modpara.def", "Lanczos_max", 18)
        rc.run(continued, "continued")
        actual = rc.rows(continued, "zvo_DynamicalGreen.dat")
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)
        tm = (continued / "output/zvo_TMComponents.dat").read_text().split()
        assert int(tm[0]) == 20, (mode, tm[:3])


cg_legacy_additional_iterations()
for tpq_method in ("TPQ", "cTPQ"):
    tpq_legacy_compat(tpq_method)
corrupt_v3_is_authoritative()
spectrum_legacy_restart()
print("restart legacy compatibility: ok (np={})".format(rc.mpi_size()),
      flush=True)

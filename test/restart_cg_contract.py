"""CG continuation reproduces the uninterrupted iteration history; periodic/converged/upper-limit checkpoints behave as specified."""
import re
import struct
import numpy as np
import restart_common as rc
ROOT = rc.fresh_dir("restart_cg_contract")
STAN = 'L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\nexct=2\nLanczos_max=40\nLanczosEps=10\noutputmode="None"\n'

def prepare(name, restart, interval=0, lanczos_max=40, eps=10):
    path = ROOT / name; rc.generate(path, STAN)
    rc.set_key(path, "calcmod.def", "ReStart", restart)
    rc.set_key(path, "modpara.def", "Lanczos_max", lanczos_max)
    rc.set_key(path, "modpara.def", "LanczosEps", eps)
    if interval:
        rc.set_key(path, "modpara.def", "CheckpointInterval", interval)
    return path

def steps(path):
    return np.loadtxt(path / "output/zvo_Lanczos_Step.dat", skiprows=1, ndmin=2)

def energies(path):
    return np.array([float(l.split()[1]) for l in (path / "output/zvo_energy.dat").read_text().splitlines() if l.lstrip().startswith("Energy")])

def same_history(a, b):
    assert a.shape == b.shape, (a.shape, b.shape)
    np.testing.assert_array_equal(a[:, 0], b[:, 0])
    np.testing.assert_allclose(a[:, 3:], b[:, 3:], rtol=0, atol=1e-9)     # Ritz values
    np.testing.assert_allclose(a[:, 1], b[:, 1], rtol=1e-6, atol=1e-12)   # residual norms

def timer_value(path, timer_id, label):
    pattern = r"^.*\[{:04d}\]\s+([0-9.eE+-]+)$".format(timer_id)
    matches = [re.match(pattern, line) for line in
               (path / "output/CalcTimer.dat").read_text().splitlines()]
    values = [float(match.group(1)) for match in matches if match]
    assert len(values) == 1 and values[0] > 0.0, (label, timer_id, values)
    print("{} timer [{}] {:.5f}".format(label, timer_id, values[0]), flush=True)
    return values[0]

# (A) uninterrupted: converges at row K (< 40) for this fixture
ref = prepare("ref", 0); rc.run(ref, "ref"); ref_steps = steps(ref); ref_energy = energies(ref)
K = int(ref_steps[-1, 0]); assert K < 40, K
# (B) periodic checkpoint every 5 iterations, stop right after the checkpoint at k=10 (STEP=11, FLAGS=0)
cont = prepare("cont", 1, interval=5)
rc.run(cont, "interrupted", binary=rc.probe(), env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:11"})
meta, rows = rc.manifest(cont)
assert meta["method"] == "CG" and rows[0]["step"] == 11 and rows[1]["step"] == 11 and rows[0]["ncomp"] == 2, rows
for k in (0, 1):
    for comp in (0, 1):
        files = rc.restart_files(cont, "zvo", k, comp, rows[k]["gen"]); assert files
        words, payload = rc.v3_header(files[0])
        assert words[22] == 3 and words[23] == k and words[24] == 11 and words[30] == 0 and words[31] == comp and words[35] == 2 and words[39] == 0
assert steps(cont)[-1, 0] == 10
# (C) continue to the end
rc.set_key(cont, "calcmod.def", "ReStart", 3); rc.drop_key(cont, "modpara.def", "CheckpointInterval")
rc.run(cont, "continued"); same_history(steps(cont), ref_steps)
np.testing.assert_allclose(energies(cont), ref_energy, rtol=0, atol=1e-10)
timer_value(cont, 7110, "restart read")
timer_value(cont, 7120, "restart output recovery")
# (D) periodic checkpoint just before convergence: stop after k=K-1 (STEP=K); resume prints row K and converges
edge = prepare("edge", 1, interval=1)
rc.run(edge, "interrupted", binary=rc.probe(), env={"HPHI_TEST_RESTART_ABORT_AFTER": "0:{}".format(K)})
rc.set_key(edge, "calcmod.def", "ReStart", 3); rc.drop_key(edge, "modpara.def", "CheckpointInterval")
rc.run(edge, "resumed"); same_history(steps(edge), ref_steps)
# (E) converged terminal checkpoint (STEP=K, FLAGS=1), same eps: no new rows, same energies, manifest untouched
conv = prepare("conv", 1); rc.run(conv, "first")
timer_value(conv, 7100, "restart write")
assert rc.manifest(conv)[1][0]["step"] == K
words, _ = rc.v3_header(rc.restart_files(conv, "zvo", 0, 0, rc.manifest(conv)[1][0]["gen"])[0]); assert words[39] == 1
before = (conv / "output/zvo_Lanczos_Step.dat").read_bytes(); manifest_before = (conv / "output/zvo_restart_manifest.dat").read_bytes()
rc.set_key(conv, "calcmod.def", "ReStart", 3); text = rc.run(conv, "noop")
assert "no remaining steps" in text and (conv / "output/zvo_Lanczos_Step.dat").read_bytes() == before
assert (conv / "output/zvo_restart_manifest.dat").read_bytes() == manifest_before
np.testing.assert_allclose(energies(conv), ref_energy, rtol=0, atol=1e-10)
# (F) upper-limit checkpoint (Lanczos_max=5 -> STEP=6, FLAGS=0) continues at iteration 6 without skipping
lim = prepare("lim", 1, lanczos_max=5); rc.run(lim, "first")
assert rc.manifest(lim)[1][0]["step"] == 6
rc.set_key(lim, "calcmod.def", "ReStart", 3); rc.set_key(lim, "modpara.def", "Lanczos_max", 40); rc.run(lim, "extended")
same_history(steps(lim), ref_steps)
# (G) upper-limit no-op in both modes: re-evaluate saved vectors, preserve STEP/history/checkpoint
for mode in (2, 3):
    lim2 = prepare("lim2_{}".format(mode), 1, lanczos_max=5); rc.run(lim2, "first")
    history = (lim2 / "output/zvo_Lanczos_Step.dat").read_bytes()
    checkpoint = {p.name: p.read_bytes() for p in (lim2 / "output").glob("zvo_restart_*")}
    (lim2 / "output/zvo_energy.dat").write_text("stale sentinel\n")
    rc.set_key(lim2, "calcmod.def", "ReStart", mode); rc.set_key(lim2, "calcmod.def", "OutputEigenVec", 1)
    text = rc.run(lim2, "noop"); assert "no remaining steps" in text
    assert np.isfinite(energies(lim2)).all() and len(energies(lim2)) == 2
    assert (lim2 / "output/zvo_Lanczos_Step.dat").read_bytes() == history
    assert {p.name: p.read_bytes() for p in (lim2 / "output").glob("zvo_restart_*")} == checkpoint
    for p in (lim2 / "output").glob("zvo_eigenvec_*_rank_*.dat"):
        assert struct.unpack_from("=i", p.read_bytes())[0] == 6
    assert list((lim2 / "output").glob("zvo_eigenvec_*_rank_*.dat"))
    # Evaluate those same output vectors through the existing import path as an independent caller check.
    evaluated = energies(lim2).copy()
    rc.set_key(lim2, "calcmod.def", "ReStart", 0); rc.set_key(lim2, "calcmod.def", "OutputEigenVec", 0)
    rc.set_key(lim2, "calcmod.def", "InputEigenVec", 1); rc.run(lim2, "import_eval")
    np.testing.assert_allclose(energies(lim2), evaluated, rtol=0, atol=1e-10)
# (H) tighter eps after a converged checkpoint: resumes and the final energy agrees with an uninterrupted tight run (plan proposal 1)
tight_ref = prepare("tight_ref", 0, eps=14); rc.run(tight_ref, "ref")
tight = prepare("tight", 1); rc.run(tight, "loose")
rc.set_key(tight, "calcmod.def", "ReStart", 3); rc.set_key(tight, "modpara.def", "LanczosEps", 14); text = rc.run(tight, "tight")
assert "no remaining steps" not in text and steps(tight)[-1, 0] > K
np.testing.assert_allclose(energies(tight), energies(tight_ref), rtol=0, atol=1e-9)
# (I) negatives: partial unit; P-only SCALAR_A tamper; mixed gen across states; missing manifest with ReStart=3; ReStart=2 fresh
broken = prepare("broken", 1); rc.run(broken, "first")
rows = rc.manifest(broken)[1]; rc.restart_files(broken, "zvo", 1, 1, rows[1]["gen"])[0].unlink()
rc.set_key(broken, "calcmod.def", "ReStart", 3); rc.run(broken, "partial", fail="restart checkpoint header scan failed")
tamper = prepare("tamper", 1); rc.run(tamper, "first"); rows = rc.manifest(tamper)[1]
for p in rc.restart_files(tamper, "zvo", 0, 1, rows[0]["gen"]):       # P files: change SCALAR_A (word 36) only
    data = bytearray(p.read_bytes()); struct.pack_into("<d", data, 8*36, 123.0); p.write_bytes(data)
rc.set_key(tamper, "calcmod.def", "ReStart", 3); rc.run(tamper, "scalar", fail="component scalars differ")
mixed = prepare("mixed", 1); rc.run(mixed, "first"); rows = rc.manifest(mixed)[1]
m = mixed / "output/zvo_restart_manifest.dat"
lines = m.read_text().splitlines()
lines = [l.replace(" gen {} ".format(rows[1]["gen"]), " gen {} ".format(1 - rows[1]["gen"])) if l.startswith("set 1 ") else l for l in lines]
m.write_text("\n".join(lines) + "\n")
rc.set_key(mixed, "calcmod.def", "ReStart", 3); rc.run(mixed, "mixed", fail="CG restart unit is incomplete or mixed")
empty = prepare("empty", 3); rc.run(empty, "missing", fail="restart manifest is required for ReStart=3")
empty2 = prepare("empty2", 2); text = rc.run(empty2, "fresh"); assert "Start from scratch" in text
# Spinless Ncond vs Nup/Ndown are the same basis (V3). Add to restart_cg_contract.py.
spinless = ROOT / "spinless_notation"
rc.generate(spinless, 'L=8\nmodel="SpinlessFermion"\nmethod="CG"\nlattice="chain"\nt=1\nnelec=4\nexct=1\nLanczos_max=40\noutputmode="None"\n')
rc.set_key(spinless, "calcmod.def", "ReStart", 1); rc.run(spinless, "ncond")
row = rc.manifest(spinless)[1][0]
words, _ = rc.v3_header(rc.restart_files(spinless, "zvo", 0, 0, row["gen"])[0])
assert words[29] == 1 and words[6] == words[7] == words[32] == 0
energy_before = energies(spinless).copy()
rc.drop_key(spinless, "modpara.def", "Ncond")
rc.set_key(spinless, "modpara.def", "Nup", 4); rc.set_key(spinless, "modpara.def", "Ndown", 0)
rc.set_key(spinless, "calcmod.def", "ReStart", 2); rc.run(spinless, "nup_ndown")
np.testing.assert_allclose(energies(spinless), energy_before, rtol=0, atol=1e-10)
# Publish from the alternate spelling, then restore Ncond and read it back.
rc.set_key(spinless, "calcmod.def", "ReStart", 1); rc.run(spinless, "publish_nup")
rc.drop_key(spinless, "modpara.def", "Nup"); rc.drop_key(spinless, "modpara.def", "Ndown")
rc.set_key(spinless, "modpara.def", "Ncond", 4); rc.set_key(spinless, "calcmod.def", "ReStart", 3)
rc.run(spinless, "restore_ncond")
np.testing.assert_allclose(energies(spinless), energy_before, rtol=0, atol=1e-10)
print('raw brief A-I and spinless notation passed: stopped row 10, resumed rows 11..{}'.format(K), flush=True)

# Additional all-state classification, reader mutations and fixed-layout sectors.
import shutil

def clone(source, name):
    path = ROOT / name
    shutil.copytree(source, path)
    return path


def checkpoint_bytes(path):
    return {p.name: p.read_bytes() for p in (path / 'output').glob('zvo_restart_*')}


def reject_preserving_history(path, label, fail='Error:', env=None):
    before = (path / 'output/zvo_Lanczos_Step.dat').read_bytes()
    saved = checkpoint_bytes(path)
    rc.run(path, label, fail=fail, env=env)
    assert (path / 'output/zvo_Lanczos_Step.dat').read_bytes() == before, label
    assert checkpoint_bytes(path) == saved, label


def change_word(path, word, value, state=0, comp=None):
    row = rc.manifest(path)[1][state]
    for component in ((0, 1) if comp is None else (comp,)):
        for f in rc.restart_files(path, 'zvo', state, component, row['gen']):
            data = bytearray(f.read_bytes())
            struct.pack_into('<Q', data, word*8, value)
            f.write_bytes(data)


# A resumed ReStart=2 run publishes at the next periodic boundary (STEP=16).
again = clone(cont, 'continue_inout')
rc.set_key(again, 'calcmod.def', 'ReStart', 2)
rc.set_key(again, 'modpara.def', 'CheckpointInterval', 5)
rc.run(again, 'second_interruption', binary=rc.probe(), env={'HPHI_TEST_RESTART_ABORT_AFTER': '0:16'})
assert rc.manifest(again)[1][0]['step'] == 16 and steps(again)[-1, 0] == 15
rc.run(again, 'finished')
same_history(steps(again), ref_steps)
np.testing.assert_allclose(energies(again), ref_energy, rtol=0, atol=1e-10)

for mode in (2, 3):
    for missing in (0, 1):
        partial = clone(conv, 'missing_state{}_mode{}'.format(missing, mode))
        m = partial / 'output/zvo_restart_manifest.dat'
        m.write_text(''.join(l for l in m.read_text().splitlines(True) if not l.startswith('set {} '.format(missing))))
        rc.set_key(partial, 'calcmod.def', 'ReStart', mode)
        reject_preserving_history(partial, 'partial', 'CG restart unit is incomplete or mixed')

# Construct the existing native 12-byte-header legacy files from valid v3 X.
legacy = clone(conv, 'legacy_complete')
for state, row in rc.manifest(legacy)[1].items():
    for f in rc.restart_files(legacy, 'zvo', state, 0, row['gen']):
        words, vector = rc.v3_header(f)
        raw = struct.pack('=iQ', words[24], len(vector)) + np.r_[0j, vector].astype('=c16').tobytes()
        (legacy / 'output/tmpvec_set{}_rank_{}.dat'.format(state, words[12])).write_bytes(raw)
(legacy / 'output/zvo_restart_manifest.dat').unlink()
rc.set_key(legacy, 'calcmod.def', 'ReStart', 3)
rc.run(legacy, 'legacy')
assert steps(legacy)[0, 0] == 1
np.testing.assert_allclose(energies(legacy), ref_energy, rtol=0, atol=1e-9)
for missing in (0, 1):
    partial = clone(legacy, 'legacy_missing{}'.format(missing))
    for f in (partial / 'output').glob('tmpvec_set{}_rank_*.dat'.format(missing)):
        f.unlink()
    reject_preserving_history(partial, 'partial', 'CG restart unit is incomplete or mixed')

mixed_legacy = clone(conv, 'mixed_legacy')
for f in (legacy / 'output').glob('tmpvec_set1_rank_*.dat'):
    shutil.copyfile(f, mixed_legacy / 'output' / f.name)
m = mixed_legacy / 'output/zvo_restart_manifest.dat'
m.write_text(''.join(l for l in m.read_text().splitlines(True) if not l.startswith('set 1 ')))
reject_preserving_history(mixed_legacy, 'mixed', 'CG restart unit is incomplete or mixed')

mixed_ranks = clone(conv, 'mixed_ranks')
m = mixed_ranks / 'output/zvo_restart_manifest.dat'
old_ranks = rc.manifest(mixed_ranks)[1][1]['ranks']
m.write_text(''.join(l.replace('ranks {} '.format(old_ranks), 'ranks {} '.format(old_ranks+1))
                     if l.startswith('set 1 ') else l for l in m.read_text().splitlines(True)))
reject_preserving_history(mixed_ranks, 'mixed', 'CG restart unit is incomplete or mixed')

# Raw Spin identity spelling and same dimension/different Sz: explicitly serial.
spin = ROOT / 'spin_notation'
rc.generate(spin, STAN.replace('L=8', 'L=4').replace('2Sz=0', '2Sz=2').replace('exct=2', 'exct=1'))
rc.set_key(spin, 'calcmod.def', 'ReStart', 1)
rc.run(spin, 'sz', mpi=False)
spin_energy = energies(spin)
rc.drop_key(spin, 'modpara.def', '2Sz')
rc.set_key(spin, 'modpara.def', 'Nup', 3)
rc.set_key(spin, 'modpara.def', 'Ndown', 1)
rc.set_key(spin, 'calcmod.def', 'ReStart', 3)
rc.run(spin, 'nup_ndown', mpi=False)
np.testing.assert_allclose(energies(spin), spin_energy, rtol=0, atol=1e-10)
rc.set_key(spin, 'calcmod.def', 'ReStart', 1)
rc.run(spin, 'publish_nup', mpi=False)
rc.drop_key(spin, 'modpara.def', 'Nup'); rc.drop_key(spin, 'modpara.def', 'Ndown')
rc.set_key(spin, 'modpara.def', '2Sz', 2); rc.set_key(spin, 'calcmod.def', 'ReStart', 3)
rc.run(spin, 'restore_sz', mpi=False)
np.testing.assert_allclose(energies(spin), spin_energy, rtol=0, atol=1e-10)
rc.set_key(spin, 'modpara.def', '2Sz', -2)
before = (spin / 'output/zvo_Lanczos_Step.dat').read_bytes()
rc.run(spin, 'other_sz', mpi=False, fail='identity / unit validation')
assert (spin / 'output/zvo_Lanczos_Step.dat').read_bytes() == before


def digest_payload(files):
    """Independent FNV framing to isolate norm/finite checks from checksum errors."""
    mask = (1 << 64) - 1
    def hash_bytes(data, h=14695981039346656037):
        for b in data:
            h = ((h ^ b) * 1099511628211) & mask
        return h
    tag = hash_bytes(b'hphi-ckpt-data-v3')
    xor = total = 0
    for f in files:
        words, vector = rc.v3_header(f)
        for i, value in enumerate(vector):
            h = hash_bytes(struct.pack('<Qdd', words[14]+i+1, value.real, value.imag), tag)
            xor ^= h; total = (total+h) & mask
    for f in files:
        data = bytearray(f.read_bytes())
        struct.pack_into('<QQ', data, 8*26, xor, total)
        f.write_bytes(data)


def mutation_matrix(source, tag, env=None):
    row = rc.manifest(source)[1][0]
    h, _ = rc.v3_header(rc.restart_files(source, 'zvo', 0, 0, row['gen'])[0])
    # Apply identity edits consistently to ranks/components, reaching identity checks.
    for name, word, value in [
        ('magic', 0, 0), ('version', 1, 99), ('phase', 2, h[2]+1),
        ('encoding', 3, 64), ('model', 4, h[4]+1), ('nsite', 5, h[5]+1),
        ('nup', 6, h[6]+1), ('ndown', 7, h[7]+1), ('ne', 8, h[8]+1),
        ('full_dim', 9, h[9]+1), ('dimension', 10, h[10]+1),
        ('group', 16, h[16]^1), ('sector_count', 17, h[17]+1),
        ('sector_xor', 18, h[18]^1), ('sector_sum', 19, h[19]+1),
        ('order', 20, h[20]^1), ('local_mask', 28, h[28]^1),
        ('fixed_flags', 29, h[29]^1), ('basis_kind', 30, 1-h[30]),
        ('total_sz', 32, h[32]^2), ('site_digest', 33, h[33]^1),
        ('order_sum', 34, h[34]+1), ('nset', 35, h[35]+1),
        ('method', 22, 1), ('state', 23, 1), ('layout', 13, 99),
        ('undefined_flags', 39, 2), ('nonfinite_scalar', 36, 0x7ff0000000000000),
    ]:
        path = clone(source, tag+'_'+name)
        change_word(path, word, value)
        reject_preserving_history(path, 'reject', env=env)
    for name, word, value in [('component', 31, 0), ('component_step', 24, h[24]+1)]:
        path = clone(source, tag+'_'+name)
        change_word(path, word, value, comp=1)
        reject_preserving_history(path, 'reject', env=env)
    for name, word, value in [('state_flags', 39, 0), ('state_time', 25, 0x3ff0000000000000),
                              ('state_scalar_b', 37, 0x3ff0000000000000),
                              ('state_scalar_c', 38, 0x3ff0000000000000)]:
        path = clone(source, tag+'_'+name)
        change_word(path, word, value, state=1)
        reject_preserving_history(path, 'reject', 'different steps, flags or scalars', env)
    for name, word in [('rank', 12), ('ranks', 11), ('offset', 14), ('local_dim', 15)]:
        path = clone(source, tag+'_'+name)
        f = rc.restart_files(path, 'zvo', 0, 0, row['gen'])[-1]
        words, _ = rc.v3_header(f)
        data = bytearray(f.read_bytes()); struct.pack_into('<Q', data, word*8, words[word]+1); f.write_bytes(data)
        reject_preserving_history(path, 'reject', env=env)
    for name in ('short_header', 'short_payload', 'extra_payload', 'nonfinite', 'norm', 'digest'):
        path = clone(source, tag+'_'+name)
        files = rc.restart_files(path, 'zvo', 0, 0, row['gen'])
        f = next(p for p in files if rc.v3_header(p)[0][15] > 0)
        data = bytearray(f.read_bytes())
        if name == 'short_header': data = data[:100]
        elif name == 'short_payload': data = data[:-1]
        elif name == 'extra_payload': data += b'x'
        else:
            vector = np.frombuffer(data, dtype='<c16', offset=320)
            if name == 'nonfinite': vector[0] = complex(float('nan'), 0)
            elif name == 'norm': vector *= 2
            else: vector *= -1  # preserves norm exactly, breaks only DATA digests
        f.write_bytes(data)
        if name in ('nonfinite', 'norm'): digest_payload(files)
        reject_preserving_history(path, 'reject', env=env)
    if h[30] == 1:
        path = clone(source, tag+'_mixed_state_layout')
        change_word(path, 13, 1-h[13], state=1)
        reject_preserving_history(path, 'reject', 'CG state identities differ', env)
    # Legal provenance-only edits are accepted, with observable results preserved.
    path = clone(source, tag+'_provenance')
    for state in (0, 1):
        change_word(path, 21, h[21]^1, state=state)
    text = rc.run(path, 'accepted', env=env)
    assert 'hamiltonian_changed='+('unknown' if h[30] == 0 else 'yes') in text
    np.testing.assert_allclose(energies(path), energies(source), rtol=0, atol=1e-10)
    print('mutation matrix {}: identity/unit/rank/length/finite/norm/digest rejected; provenance accepted'.format(tag), flush=True)

mutation_matrix(conv, 'raw')

import symmetry_general_terms as sg
sg.ROOT = ROOT / 'sector_fixtures'
sg.ROOT.mkdir()


def sector_case(template, model, length, momentum, states, raw, projector):
    rc.set_key(template, 'mod.def', 'exct', 2)
    rc.set_key(template, 'mod.def', 'Lanczos_max', 100)
    rc.set_key(template, 'mod.def', 'LanczosEps', 10)
    (template / 'namelist.def').write_bytes((template / 'sym.def').read_bytes())
    weights, vectors = np.linalg.eigh(projector)
    basis = vectors[:, weights > .5]
    exact = np.linalg.eigvalsh(basis.conj().T @ raw @ basis)[:2]
    for layout in ('replicated', 'distributed'):
        prefix = 'sector_'+layout
        env = {'HPHI_SYMMETRY_BASIS_LAYOUT': layout}
        ref = clone(template, prefix+'_ref')
        rc.set_key(ref, 'calc.def', 'ReStart', 0)
        rc.run(ref, 'uninterrupted', env=env)
        history = steps(ref)
        K = int(history[-1, 0]); assert K > 3
        np.testing.assert_allclose(energies(ref), exact, rtol=0, atol=1e-8)
        cont = clone(template, prefix+'_cont')
        rc.set_key(cont, 'calc.def', 'ReStart', 1)
        rc.set_key(cont, 'mod.def', 'CheckpointInterval', 2)
        rc.run(cont, 'interrupted', binary=rc.probe(), env=dict(env, HPHI_TEST_RESTART_ABORT_AFTER='0:3'))
        assert steps(cont)[-1, 0] == 2 and rc.manifest(cont)[1][0]['step'] == 3
        rc.set_key(cont, 'calc.def', 'ReStart', 2)
        text = rc.run(cont, 'continued', env=env)
        assert 'no remaining steps' not in text
        same_history(steps(cont), history)
        np.testing.assert_allclose(energies(cont), energies(ref), rtol=0, atol=1e-10)
        row = rc.manifest(cont)[1][0]
        w, _ = rc.v3_header(rc.restart_files(cont, 'zvo', 0, 0, row['gen'])[0])
        assert row['step'] == K and w[39] == 1 and w[30] == 1
        # Converged no-op: remove injected tail once, then remain byte invariant.
        clean = (cont / 'output/zvo_Lanczos_Step.dat').read_bytes()
        with (cont / 'output/zvo_Lanczos_Step.dat').open('ab') as stream:
            stream.write(('{} 0 0 0 0\npartial'.format(K+1)).encode())
        saved = checkpoint_bytes(cont)
        rc.set_key(cont, 'mod.def', 'Lanczos_max', 1)
        rc.drop_key(cont, 'mod.def', 'CheckpointInterval')
        for mode in (2, 3):
            rc.set_key(cont, 'calc.def', 'ReStart', mode)
            (cont / 'output/zvo_energy.dat').write_text('stale\n')
            assert 'no remaining steps' in rc.run(cont, 'converged_noop{}'.format(mode), env=env)
            assert (cont / 'output/zvo_Lanczos_Step.dat').read_bytes() == clean
            assert checkpoint_bytes(cont) == saved
            np.testing.assert_allclose(energies(cont), energies(ref), rtol=0, atol=1e-10)
        for mode in (2, 3):
            lim = clone(template, prefix+'_limit{}'.format(mode))
            rc.set_key(lim, 'calc.def', 'ReStart', 1)
            rc.set_key(lim, 'mod.def', 'Lanczos_max', 2)
            rc.run(lim, 'limited', env=env)
            assert rc.manifest(lim)[1][0]['step'] == 3
            clean = (lim / 'output/zvo_Lanczos_Step.dat').read_bytes()
            saved = checkpoint_bytes(lim)
            with (lim / 'output/zvo_Lanczos_Step.dat').open('ab') as stream:
                stream.write(b'3 0 0 0 0\npartial')
            rc.set_key(lim, 'calc.def', 'ReStart', mode)
            rc.set_key(lim, 'calc.def', 'OutputEigenVec', 1)
            rc.set_key(lim, 'mod.def', 'Lanczos_max', 1)
            for repeat in (0, 1):
                (lim / 'output/zvo_energy.dat').write_text('stale\n')
                assert 'no remaining steps' in rc.run(lim, 'limit_noop{}'.format(repeat), env=env)
                assert (lim / 'output/zvo_Lanczos_Step.dat').read_bytes() == clean
                assert checkpoint_bytes(lim) == saved
                assert len(energies(lim)) == 2 and np.isfinite(energies(lim)).all()
                if repeat == 0: evaluated = energies(lim)
                else: np.testing.assert_array_equal(energies(lim), evaluated)
            outputs = list((lim / 'output').glob('zvo_eigenvec_*_rank_*.dat'))
            assert outputs
            for f in outputs:
                data = f.read_bytes()
                assert struct.unpack_from('<Q', data, 24*8)[0] == 3
                assert struct.unpack_from('<Q', data, 8)[0] in (1, 2)
            rc.set_key(lim, 'calc.def', 'ReStart', 0)
            rc.set_key(lim, 'calc.def', 'OutputEigenVec', 0)
            rc.set_key(lim, 'calc.def', 'InputEigenVec', 1)
            rc.run(lim, 'import_eval', env=env)
            np.testing.assert_allclose(energies(lim), evaluated, rtol=0, atol=1e-10)
        mutation_matrix(cont, prefix, env)
        dims = [rc.v3_header(f)[0][15] for f in rc.restart_files(cont, 'zvo', 0, 0, row['gen'])]
        print('{}: stopped row 2, resumed rows 3..{}, dim={}, empty owners={}'.format(prefix, K, sum(dims), dims.count(0)), flush=True)

sg.prepare('Spin', 8, sector_test=sector_case, sector_momenta=(0,))
print('CG raw/sector contract and reader matrix ok', flush=True)

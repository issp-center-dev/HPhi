"""Portable checkpoint payloads versus independent sector eigenvectors.

Every negative load corrupts just one rank and must terminate collectively.
"""
import os
from pathlib import Path
import shutil
import struct

import numpy as np
import symmetry_general_terms as fixture

fixture.ROOT = Path("symmetry_checkpoint")
fixture.ROOT.mkdir(exist_ok=True)
WORDS = 28


def read_file(path):
    data = path.read_bytes()
    header = struct.unpack("<{}Q".format(WORDS), data[:8*WORDS])
    assert header[:4] == (0x0a31565349485048, 1, 1, 128)
    assert len(data) == 8*WORDS + 16*header[15]
    return header, np.frombuffer(data[8*WORDS:], dtype="<c16")


def check_sector(path, model, length, momentum, states, raw, projector):
    columns, covered = [], set()
    for index in range(len(states)):
        if index in covered or np.linalg.norm(projector[:, index]) < 1e-12:
            continue
        col = projector[:, index]
        covered.update(np.flatnonzero(abs(col) > 1e-12))
        columns.append(col / np.linalg.norm(col))
    basis = np.column_stack(columns)
    h = basis.conj().T @ raw @ basis
    dim = h.shape[0]
    calc = (path / "calc.def").read_text()
    for layout in ("replicated", "distributed"):
        output = path / "output"
        if output.exists():
            shutil.rmtree(output)
        (path / "calc.def").write_text(calc + "OutputEigenVec 1\n")
        energy = fixture.run(path, "write_k{}_{}".format(momentum, layout), layout=layout)
        np.testing.assert_allclose(energy, np.linalg.eigvalsh(h), atol=3e-8, rtol=0)
        files = sorted(output.glob("zvo_eigenvec_0_rank_*.dat"),
                       key=lambda p: int(p.stem.split("rank_")[1]))
        assert files
        for state in range(dim):
            parts = []
            offset = 0
            for rank in range(len(files)):
                header, vector = read_file(output / "zvo_eigenvec_{}_rank_{}.dat".format(state, rank))
                assert header[9:12] == (len(states), dim, len(files))
                assert header[12:15] == (rank, int(layout == "distributed"), offset)
                assert header[22:24] == (3, state)
                assert header[25] == 0  # CG physical time
                parts.extend(vector)
                offset += header[15]
            vector = np.array(parts)
            assert offset == dim
            np.testing.assert_allclose(np.vdot(vector, vector), 1, atol=1e-10)
            np.testing.assert_allclose(h @ vector, energy[state]*vector, atol=3e-7, rtol=0)
        saved_files = {p: p.read_bytes() for p in output.glob("*eigenvec*.dat")}
        (path / "calc.def").write_text(calc + "InputEigenVec 1\nOutputEigenVec 1\n")
        np.testing.assert_allclose(fixture.run(path, "read_k{}_{}".format(momentum, layout), layout=layout),
                                   energy, atol=2e-12, rtol=0)
        for p, data in saved_files.items():
            assert p.read_bytes() == data, "checkpoint round-trip changed bytes"
        (path / "calc.def").write_text(calc + "InputEigenVec 1\n")
        # Group-row order does not change vector index or phase conventions.
        transsym = path / "group.def"
        original_group = transsym.read_text()
        lines = original_group.splitlines()
        transsym.write_text("\n".join(lines[:5] + lines[5:5+length][::-1] + lines[5+length:][::-1]) + "\n")
        np.testing.assert_allclose(fixture.run(path, "reorder_k{}_{}".format(momentum, layout), layout=layout),
                                   energy, atol=2e-12, rtol=0)
        transsym.write_text(original_group)
        shutil.copytree(output, path / "k{}_{}".format(momentum, layout))
        if model != "Spin" or momentum != 0:
            continue
        (path / "calc.def").write_text(calc + "InputEigenVec 2\n")
        fixture.run(path, "reject_text_input_{}".format(layout), layout=layout,
                    fail="InputEigenVec=1 only")
        (path / "calc.def").write_text(calc + "InputEigenVec 1\n")
        # k=0 and k=pi have the same dimension for this Spin L6 sector.
        # A real different character must still be rejected, not just bad sizes.
        lines = original_group.splitlines()
        for row in range(5, 5+length):
            group_id = int(lines[row].split()[0])
            lines[row] = "{} {} 0".format(group_id, (-1)**group_id)
        transsym.write_text("\n".join(lines) + "\n")
        fixture.run(path, "reject_other_sector_{}".format(layout), layout=layout,
                    fail="checkpoint header / sector / layout validation failed")
        transsym.write_text(original_group)
        original = {p: p.read_bytes() for p in output.glob("*eigenvec*.dat")}
        victim = files[-1]
        saved = victim.read_bytes()
        for word in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20):
            damaged = bytearray(saved)
            damaged[8*word] ^= 1
            victim.write_bytes(damaged)
            fixture.run(path, "reject_header_{}_{}".format(layout, word), layout=layout,
                        fail="checkpoint header / sector / layout validation failed")
            victim.write_bytes(saved)
        for label, payload in [("truncated_header", saved[:32]),
                               ("truncated_payload", saved[:-1]),
                               ("trailing", saved + b"x")]:
            victim.write_bytes(payload)
            fixture.run(path, "reject_{}_{}".format(layout, label), layout=layout,
                        fail="checkpoint")
            victim.write_bytes(saved)
        victim.unlink()
        fixture.run(path, "reject_missing_{}".format(layout), layout=layout,
                    fail="checkpoint header / sector / layout validation failed")
        victim.write_bytes(saved)
        damaged = bytearray(saved)
        damaged[8*23] ^= 1  # different state index, even with identical dimension
        victim.write_bytes(damaged)
        fixture.run(path, "reject_state_{}".format(layout), layout=layout, fail="checkpoint")
        victim.write_bytes(saved)
        owned = files[0]
        data = bytearray(owned.read_bytes())
        data[8*WORDS] ^= 1
        owned.write_bytes(data)
        fixture.run(path, "reject_checksum_{}".format(layout), layout=layout,
                    fail="checkpoint payload checksum failed")
        owned.write_bytes(original[owned])
        data = bytearray(owned.read_bytes())
        data[8*WORDS:8*WORDS+8] = struct.pack("<d", float("nan"))
        owned.write_bytes(data)
        fixture.run(path, "reject_nonfinite_{}".format(layout), layout=layout,
                    fail="checkpoint finite vector validation failed")
        owned.write_bytes(original[owned])
        for p in files:
            data = p.read_bytes()
            p.write_bytes(data[:8*WORDS] + bytes(len(data)-8*WORDS))
        fixture.run(path, "reject_zero_norm_{}".format(layout), layout=layout,
                    fail="checkpoint global norm validation failed")
        for p in files:
            p.write_bytes(original[p])
        fixture.run(path, "reject_layout_{}".format(layout),
                    layout="distributed" if layout == "replicated" else "replicated",
                    fail="checkpoint header / sector / layout validation failed")
        # A different Hamiltonian is valid for a quench-capable state import.
        # A uniform fixed-Sz field changes H by a known scalar.
        transfer = path / "Trans.def"
        before = transfer.read_text()
        rows = before.splitlines()
        for i in range(5, len(rows)):
            fields = rows[i].split()
            fields[4] = str(float(fields[4]) + .2)
            rows[i] = " ".join(fields)
        transfer.write_text("\n".join(rows) + "\n")
        # The fixture's Spin diagonal Trans is n_down, sum fixed at L/2.
        np.testing.assert_allclose(fixture.run(path, "quench_{}".format(layout), layout=layout),
                                   energy - .2*(length//2), atol=2e-12, rtol=0)
        assert "hamiltonian_changed=yes" in (path / "quench_{}.log".format(layout)).read_text()
        transfer.write_text(before)
        # A failed final rename on one rank must be reported by the CG driver.
        victim.unlink()
        victim.mkdir()
        (path / "calc.def").write_text(calc + "OutputEigenVec 1\n")
        fixture.run(path, "reject_write_{}".format(layout), layout=layout,
                    fail="checkpoint publish failed")
        victim.rmdir()
        for p, content in original.items():
            p.write_bytes(content)
        blocked = Path(str(victim) + ".part")
        blocked.mkdir()
        fixture.run(path, "reject_open_{}".format(layout), layout=layout,
                    fail="checkpoint write failed")
        assert blocked.is_dir()  # an unsuccessful fopen must not remove it
        blocked.rmdir()
    print("{} k={}: checkpoint eigenvectors, reload, and identity checks PASS".format(model, momentum))


for model, length in [("Spin", 6), ("SpinlessFermion", 6), ("Hubbard", 4), ("tJ", 4)]:
    fixture.prepare(model, length, sector_test=check_sector)

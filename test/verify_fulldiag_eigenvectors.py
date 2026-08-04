#!/usr/bin/env python3
"""Validate FullDiag eigenvector files and optionally compare eigenspaces."""

from __future__ import annotations

import glob
import math
import os
import struct
import sys
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_eigenvalues(output_dir: Path) -> list[float]:
    path = output_dir / "Eigenvalue.dat"
    values: list[float] = []
    try:
        with path.open(encoding="utf-8") as stream:
            for line in stream:
                fields = line.split()
                if len(fields) >= 2:
                    values.append(float(fields[1]))
    except OSError as exc:
        fail(f"cannot read {path}: {exc}")
    if not values:
        fail(f"{path} contains no eigenvalues")
    return values


def read_vectors(output_dir: Path, head: str) -> list[list[complex]]:
    pattern = str(output_dir / f"{head}_eigenvec_*_rank_0.dat")
    paths = glob.glob(pattern)
    indexed: dict[int, Path] = {}
    for name in paths:
        path = Path(name)
        prefix = f"{head}_eigenvec_"
        suffix = "_rank_0.dat"
        stem = path.name
        if not stem.startswith(prefix) or not stem.endswith(suffix):
            fail(f"unexpected FullDiag eigenvector name: {path}")
        try:
            state = int(stem[len(prefix) : -len(suffix)])
        except ValueError:
            fail(f"invalid state index in {path}")
        indexed[state] = path

    other_rank_files = glob.glob(
        str(output_dir / f"{head}_eigenvec_*_rank_[1-9]*.dat")
    )
    if other_rank_files:
        fail(
            "FullDiag emitted rank-split files instead of global rank_0 files: "
            + ", ".join(other_rank_files)
        )
    if sorted(indexed) != list(range(len(indexed))):
        fail(f"non-contiguous FullDiag eigenvector indices: {sorted(indexed)}")

    vectors: list[list[complex]] = []
    dimension: int | None = None
    int_size = struct.calcsize("@i")
    ulong_size = struct.calcsize("@L")
    complex_size = struct.calcsize("@dd")
    for state in range(len(indexed)):
        path = indexed[state]
        data = path.read_bytes()
        header_size = int_size + ulong_size
        if len(data) < header_size:
            fail(f"truncated header in {path}")
        iteration = struct.unpack("@i", data[:int_size])[0]
        local_dimension = struct.unpack(
            "@L", data[int_size : int_size + ulong_size]
        )[0]
        expected_size = header_size + (local_dimension + 1) * complex_size
        if len(data) != expected_size:
            fail(
                f"incorrect size for {path}: got {len(data)}, "
                f"expected {expected_size}"
            )
        if iteration != 0:
            fail(f"FullDiag iteration header in {path} is {iteration}, expected 0")
        if dimension is None:
            dimension = local_dimension
        elif local_dimension != dimension:
            fail(
                f"dimension mismatch in {path}: {local_dimension} != {dimension}"
            )
        components = [
            complex(*value)
            for value in struct.iter_unpack("@dd", data[header_size:])
        ]
        if abs(components[0]) > 1.0e-15:
            fail(f"unused first component in {path} is not zero")
        vectors.append(components[1:])

    if dimension is None or len(vectors) != dimension:
        fail(
            f"expected one vector per Hilbert-space dimension, "
            f"got {len(vectors)} vectors of dimension {dimension}"
        )
    return vectors


def check_orthonormal(vectors: list[list[complex]], tolerance: float) -> None:
    dimension = len(vectors)
    max_error = 0.0
    for i in range(dimension):
        for j in range(dimension):
            overlap = sum(
                vectors[i][k].conjugate() * vectors[j][k]
                for k in range(dimension)
            )
            target = 1.0 if i == j else 0.0
            max_error = max(max_error, abs(overlap - target))
    if not math.isfinite(max_error) or max_error > tolerance:
        fail(f"eigenvectors are not orthonormal: max error = {max_error:.3e}")


def degenerate_groups(values: list[float], tolerance: float) -> list[range]:
    groups: list[range] = []
    begin = 0
    for end in range(1, len(values) + 1):
        if end == len(values) or abs(values[end] - values[begin]) > tolerance:
            groups.append(range(begin, end))
            begin = end
    return groups


def compare_eigenspaces(
    values: list[float],
    vectors: list[list[complex]],
    reference_values: list[float],
    reference_vectors: list[list[complex]],
    tolerance: float,
) -> None:
    if len(values) != len(reference_values):
        fail("eigenvalue counts differ between FullDiag solver outputs")
    max_eigenvalue_error = max(
        abs(value - reference)
        for value, reference in zip(values, reference_values)
    )
    if max_eigenvalue_error > tolerance:
        fail(
            "eigenvalues differ between FullDiag solver outputs: "
            f"max error = {max_eigenvalue_error:.3e}"
        )

    dimension = len(values)
    max_projector_error = 0.0
    for group in degenerate_groups(reference_values, tolerance):
        for row in range(dimension):
            for column in range(dimension):
                projector = sum(
                    vectors[state][row] * vectors[state][column].conjugate()
                    for state in group
                )
                reference_projector = sum(
                    reference_vectors[state][row]
                    * reference_vectors[state][column].conjugate()
                    for state in group
                )
                max_projector_error = max(
                    max_projector_error, abs(projector - reference_projector)
                )
    if not math.isfinite(max_projector_error) or max_projector_error > tolerance:
        fail(
            "FullDiag eigenspaces differ from the LAPACK reference: "
            f"max projector error = {max_projector_error:.3e}"
        )


def load(output_dir: Path, head: str) -> tuple[list[float], list[list[complex]]]:
    values = read_eigenvalues(output_dir)
    vectors = read_vectors(output_dir, head)
    if len(values) != len(vectors):
        fail(
            f"eigenvalue/vector count mismatch in {output_dir}: "
            f"{len(values)} != {len(vectors)}"
        )
    check_orthonormal(vectors, 1.0e-10)
    return values, vectors


def main() -> None:
    if len(sys.argv) not in (3, 4):
        fail(
            "usage: verify_fulldiag_eigenvectors.py OUTPUT_DIR HEAD "
            "[REFERENCE_OUTPUT_DIR]"
        )
    output_dir = Path(sys.argv[1])
    head = sys.argv[2]
    values, vectors = load(output_dir, head)
    if len(sys.argv) == 4:
        reference_values, reference_vectors = load(Path(sys.argv[3]), head)
        compare_eigenspaces(
            values,
            vectors,
            reference_values,
            reference_vectors,
            2.0e-8,
        )
    print(
        f"verified {len(vectors)} FullDiag eigenvectors in "
        f"{os.fspath(output_dir)}"
    )


if __name__ == "__main__":
    main()

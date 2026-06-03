Restart Mechanism
=================

This section describes the restart mechanism in :math:`{\mathcal H}\Phi`.

Overview
--------

:math:`{\mathcal H}\Phi` supports checkpoint/restart functionality, allowing:

- Resumption of interrupted calculations
- Breaking large calculations into smaller jobs
- Iterative refinement of eigenvectors

Restart Files
-------------

The following files are used for restart:

**Eigenvector files**
  ``zvo_eigenvec_{rank}_set{set}.dat``
    Stores the eigenvector for each MPI rank.

**Restart vector files**
  ``restart_{rank}.dat``
    Contains the Lanczos vectors needed to resume iteration.

**Tri-diagonal matrix**
  ``zvo_Lanczos_Step.dat``
    Stores the Lanczos coefficients.

Input/Output Settings
---------------------

In the ``CalcMod`` input file:

``iInputEigenVec``
  - ``0``: Start from random or specified initial vector (default)
  - ``1``: Read eigenvector from file and recalculate
  - ``2``: Read eigenvector from file

``iOutputEigenVec``
  - ``0``: Do not output eigenvector (default)
  - ``1``: Output eigenvector to file

Usage Examples
--------------

**Initial calculation with eigenvector output**

Set ``iOutputEigenVec = 1`` in the CalcMod file.
After the calculation completes, eigenvector files will be written.

**Restart from eigenvector**

Set ``iInputEigenVec = 1`` and ``iOutputEigenVec = 1``.
The calculation will read the eigenvector from the previous run
and refine it.

**Calculate expectation values only**

Set ``iInputEigenVec = 2``.
The eigenvector will be read and used to compute observables
without further diagonalization.

File Format
-----------

Eigenvector files are stored in binary format for efficiency.
Each file contains:

1. Header information (system size, etc.)
2. Complex double array of eigenvector components

The files are platform-dependent due to binary format.
Ensure restart files are used on the same architecture.

MPI Considerations
------------------

When using MPI, each process writes and reads its own restart files.
The number of MPI processes must match between the original and restart runs.

If you need to change the process count, you must:

1. Gather all eigenvector files
2. Redistribute the data
3. Split into new per-process files

Common Issues
-------------

**File not found**
  Ensure restart files are in the working directory with correct names.

**Process count mismatch**
  The number of MPI processes must match the original calculation.

**File corruption**
  If calculations produce unexpected results after restart,
  the restart files may be corrupted. Delete them and start fresh.

Best Practices
--------------

- Regularly save eigenvector files during long calculations
- Keep backup copies of restart files
- Verify results after restart by checking conserved quantities
- Document the calculation parameters used with each restart file

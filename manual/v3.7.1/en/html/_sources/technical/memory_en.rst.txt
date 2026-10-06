Memory Management
=================

This section describes how :math:`{\mathcal H}\Phi` manages memory for large-scale calculations.

Memory Requirements
-------------------

The dominant memory usage in :math:`{\mathcal H}\Phi` comes from storing the state vectors.
For a Hilbert space of dimension :math:`D` distributed across :math:`P` MPI processes,
each process requires:

**Lanczos method**
  - Two state vectors (``v0``, ``v1``): :math:`2 \times D/P \times 16` bytes (complex double)
  - One buffer for MPI communication: :math:`D/P \times 16` bytes

**LOBPCG method**
  - Multiple state vectors for subspace iteration
  - Buffer for orthogonalization

**List arrays**
  - ``list_1``: :math:`D/P \times 8` bytes (unsigned long)
  - ``list_2_1``, ``list_2_2``: depends on system size

Memory Allocation
-----------------

:math:`{\mathcal H}\Phi` allocates memory dynamically using wrapper functions defined in
``common/setmemory.c``. These functions handle allocation for various data types:

- ``d_1d_allocate``: 1D double array
- ``cd_1d_allocate``: 1D complex double array
- ``i_2d_allocate``: 2D integer array
- ``lui_1d_allocate``: 1D unsigned long integer array

Memory is allocated at the beginning of the calculation in ``xsetmem.c`` and
freed at the end.

Estimating Memory Usage
-----------------------

Before running a calculation, you can estimate the required memory:

.. math::

   M_{\text{total}} \approx \frac{D}{P} \times (N_v \times 16 + N_l \times 8) \text{ bytes}

where:

- :math:`D` is the Hilbert space dimension
- :math:`P` is the number of MPI processes
- :math:`N_v` is the number of state vectors (typically 3-4 for Lanczos)
- :math:`N_l` is the number of list arrays (typically 1-3)

**Example**
  For a 16-site Hubbard model with 4 up and 4 down electrons (:math:`S_z = 0`):

  - Hilbert space: :math:`D = \binom{16}{4}^2 = 1,820^2 = 3,312,400`
  - With 4 MPI processes: :math:`D/P = 828,100` states per process
  - Memory per process: approximately 50 MB for state vectors

Memory Optimization Strategies
------------------------------

**Increasing MPI processes**
  Distributing the Hilbert space across more processes reduces
  per-process memory requirements.

**Using symmetries**
  Particle number conservation, :math:`S_z` conservation, and other symmetries
  reduce the Hilbert space dimension.

**Restart functionality**
  For very large calculations, the restart feature allows breaking
  the calculation into smaller chunks.

Common Memory Issues
--------------------

**Out-of-memory errors**
  If the calculation runs out of memory, try:

  1. Increase the number of MPI processes
  2. Use symmetry restrictions
  3. Request more memory per node on cluster systems

**Memory leaks**
  :math:`{\mathcal H}\Phi` carefully manages memory to avoid leaks.
  If you observe growing memory usage, please report it as a bug.

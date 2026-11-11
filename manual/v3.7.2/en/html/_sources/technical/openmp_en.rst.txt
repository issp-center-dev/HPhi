OpenMP Parallelization
======================

This section describes how :math:`{\mathcal H}\Phi` utilizes OpenMP for shared-memory parallelization.

Overview
--------

:math:`{\mathcal H}\Phi` uses OpenMP to parallelize operations over the local Hilbert space
on each MPI process. This hybrid MPI+OpenMP approach enables efficient use of modern
multi-core computing nodes.

Parallelized Operations
-----------------------

The following operations are parallelized using OpenMP:

**Hamiltonian-vector multiplication**
  The main computational kernel in iterative eigensolvers. Each thread processes
  a portion of the local state vector.

**Diagonal matrix elements**
  Evaluation of on-site interactions and chemical potentials.

**Inner products and norms**
  Vector operations required for Lanczos and LOBPCG algorithms.

**Expectation values**
  Physical observable calculations such as energy, spin correlations,
  and charge correlations.

Thread-safe Implementation
--------------------------

Care must be taken to avoid race conditions when multiple threads update
shared data structures. :math:`{\mathcal H}\Phi` uses the following strategies:

**Separate accumulation**
  Each thread accumulates results to separate memory locations,
  which are then combined in a serial section or using OpenMP reduction clauses.

**Read-only shared data**
  Hamiltonian parameters and state indexing arrays are read-only during
  parallel regions.

Setting the Number of Threads
-----------------------------

The number of OpenMP threads is controlled by the ``OMP_NUM_THREADS`` environment variable:

.. code-block:: bash

   export OMP_NUM_THREADS=4
   mpirun -np 2 ./HPhi -e namelist.def

For optimal performance, the product of MPI processes and OpenMP threads
should match the number of available CPU cores.

Performance Considerations
--------------------------

**Load balancing**
  The local Hilbert space is evenly divided among threads by default.
  For systems with symmetry restrictions, some threads may have more work.

**Memory bandwidth**
  Hamiltonian-vector multiplication is often memory-bound.
  Using fewer threads with better cache locality may sometimes be more efficient
  than maximizing thread count.

**NUMA effects**
  On NUMA systems, memory affinity can significantly impact performance.
  Consider using ``OMP_PROC_BIND`` and ``OMP_PLACES`` for thread pinning.

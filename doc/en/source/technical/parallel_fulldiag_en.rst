.. highlight:: none

.. _Sec:ParallelFullDiag:

Parallel full diagonalization
=============================

This section describes the parallel algorithms used by the full
diagonalization method (``method="FullDiag"``) and compares the speed of
the available backends. See :ref:`Subsec:calcmod` for the specification
of the input keywords (``Solver``, ``ExpecMode``, ``NGPU``).

Computation pipeline
--------------------

A full diagonalization run consists of four stages:

1. **Hamiltonian generation**: the :math:`N \times N` Hermitian matrix
   :math:`H_{ij} = \langle \psi_i | \hat{H} | \psi_j \rangle` is built
   in the real-space configuration basis.
2. **Diagonalization**: all eigenvalues and eigenvectors are computed by
   the backend selected with the ``Solver`` keyword (LAPACK / ScaLAPACK
   / MAGMA / ELPA).
3. **Eigenvector redistribution** (distributed ``Solver 1``/``3`` runs
   with ``ExpecMode 1``/``2``): the eigenvectors are rearranged from
   the distribution used by the solver into a distribution suited for
   observable evaluation.
4. **Observable evaluation and output**: for each eigenstate
   :math:`|\Phi_i\rangle`, expectation values
   :math:`\langle \Phi_i | \hat{A} | \Phi_i \rangle` (energy,
   fluctuations, Green functions) are computed and written. The
   evaluation kernel is selected with ``ExpecMode``.

Distributed Hamiltonian generation (Solver 3)
---------------------------------------------

With ``Solver 0``/``1``/``2`` every rank holds a replicated copy of the
full :math:`N \times N` matrix, so the per-process memory is
:math:`O(N^2)`. When ``Solver 3`` (ELPA) runs with more than one MPI
process, each rank generates and stores only its own contiguous column
block (panel) of the matrix. The peak memory per rank becomes
:math:`O(N^2/P)` (:math:`P`: number of processes), so larger dimensions
become reachable by adding processes. As a measured example, for the
8-site Hubbard chain (:math:`N=4900`) the replicated layout takes
1.93 GB in one process while the distributed layout takes 0.63 GB per
rank with four processes.

The generated panels are then redistributed by MPI communication into
the two-dimensional block-cyclic layout required by ELPA (and by the
ScaLAPACK descriptors). The process grid :math:`n_{\rm prow} \times
n_{\rm pcol}` is chosen as close to square as :math:`P` allows. If the
matrix dimension :math:`N` is smaller than
:math:`\max(n_{\rm prow}, n_{\rm pcol})`, the run stops before the ELPA
diagonalization with an error asking to reduce the number of ranks.

Diagonalization with ELPA
-------------------------

ELPA (Eigenvalue soLvers for Petaflop Applications) is a parallel
library for dense eigenvalue problems that uses the same
two-dimensional block-cyclic distribution as ScaLAPACK while achieving
higher parallel efficiency. :math:`{\mathcal H}\Phi` uses the complex
Hermitian solver and selects the 1-stage/2-stage algorithm itself: GPU
runs always use 1-stage, and CPU runs use 2-stage when the block size
can take its default value, falling back to 1-stage when the matrix is
small relative to the process grid and the block size is capped.
Setting ``NGPU`` :math:`\geq 1` enables the GPU version; the
process-to-GPU assignment is handled by ELPA (one process per GPU is
the recommended launch configuration; requires ELPA 2023.11.001 or
later). ``NGPU`` is an enablement flag and does not physically limit
the devices used — control the assignment through the job scheduler or
``CUDA_VISIBLE_DEVICES``.

State-task-parallel observable evaluation (ExpecMode 1)
-------------------------------------------------------

In the conventional evaluation (``ExpecMode 0``) the eigenstates are
processed one by one with all ranks cooperating, which requires MPI
collectives (eigenvector gather and reductions) for every state. Since
the number of states equals :math:`N`, this communication dominates the
run time, especially on multiple nodes.

With ``ExpecMode 1`` the eigenvectors are redistributed after the
diagonalization into a "state panel" layout in which each rank owns a
contiguous block of eigenstates, and each rank then evaluates its own
states independently, **with no communication inside the per-state
loop**. The physical quantities (``zvo_phys*.dat``) are collected from
the ranks only once at the end. The aggregate Green-function output
(``OutputGreenFormat 1``) is written to per-rank partial files, and
rank 0 merges them only after the manifests of all ranks report
success, publishing the final file transactionally via a temporary
file and rename (a design that avoids publishing an incomplete final
file on failure).

Trace-kernel evaluation (ExpecMode 2)
-------------------------------------

Even with ``ExpecMode 1``, the one-body and two-body Green functions
are evaluated by a loop that applies each operator state by state,
paying the per-operator basis-scan overhead :math:`N` times. With
``ExpecMode 2``, for each operator :math:`\hat{A}` the basis mapping
(target state index :math:`k \to k'` and amplitude) is precomputed
**once**, and all owned eigenstates are streamed through a dense loop
over that mapping. The per-operator overhead is thus amortized over the
number of states, which pays off when one- and two-body Green functions
are evaluated for many eigenstates.

The supported models in this version are ``Hubbard``/``HubbardGC`` and
spin-1/2 ``Spin``/``SpinGC``. For unsupported models, and under
certain runtime conditions (evaluator sharing with the N-body Green
functions, no operators of a kind defined, or the result buffer
exceeding ``HPHI_TRACE_BUF_MAX_MB``), the affected quantity
automatically falls back to the ``ExpecMode 1`` path, and the decision
is reported by ``INFO`` lines just before the observable evaluation. See the
``ExpecMode`` entry in :ref:`Subsec:calcmod` for details.

Speed comparison
----------------

Single node
~~~~~~~~~~~

A spin-1/2 chain in a transverse field (``model="SpinGC"``,
``lattice="chain"``, :math:`J=1`, :math:`\Gamma=0.5`,
:math:`L=8,10,12,14`, dimension :math:`N=2^L=256,\dots,16384`) was
measured on one node of the ISSP supercomputer kugui (AMD EPYC 7763,
Intel compilers + Intel MPI + MKL). Every CPU configuration uses 16
cores in total: ``Solver 0`` uses 1 process :math:`\times` 16 threads,
``Solver 1``/``3`` use 16 processes :math:`\times` 1 thread. The GPU
runs use a GPU node of the same system (EPYC 7763 + 2 :math:`\times`
NVIDIA A100-SXM4-40GB), with one process per GPU (``NGPU`` 1 and 2,
assigned automatically by ELPA).

Diagonalization time (``LapackDiag`` section of ``CalcTimer.dat``,
seconds):

.. csv-table::
   :header: ":math:`N`", "Solver 0 (LAPACK)", "Solver 1 (ScaLAPACK)", "Solver 3 (ELPA, CPU)", "Solver 3 (ELPA, GPU :math:`\times` 1)", "Solver 3 (ELPA, GPU :math:`\times` 2)"
   :widths: 10, 18, 18, 18, 18, 18

   "256", "0.82", "0.19", "0.20", "3.42", "1.69"
   "1024", "0.94", "0.33", "0.21", "1.40", "1.51"
   "4096", "42.22", "7.60", "3.72", "3.41", "3.11"
   "16384", "1929.6", "536.9", "172.6", "48.23", "35.51"

.. figure:: ../../../figs/fulldiag_solver_bench.png
   :name: fig_fulldiag_solver_bench_en
   :alt: Diagonalization time versus matrix dimension
   :width: 600px

   Diagonalization time versus matrix dimension :math:`N` (log-log).

At :math:`N=16384`, ``Solver 3`` (ELPA CPU, 16 processes) is about 11
times faster than ``Solver 0`` (LAPACK) and about 3.1 times faster
than ``Solver 1`` (ScaLAPACK). The GPU runs pay an initialization and
transfer overhead at small :math:`N`, break even with the 16-core CPU
runs around :math:`N \sim 4000`, and reach about 3.6 times the ELPA
CPU speed (about 40 times LAPACK) with one A100 and about 4.9 times
(about 54 times LAPACK) with two A100s at :math:`N=16384`.

Observable-evaluation time (``Solver 3`` CPU with 16 processes,
aggregate one- and two-body Green-function output for all states,
``CalcPhys`` section of ``CalcTimer.dat``, seconds):

.. csv-table::
   :header: ":math:`N`", "ExpecMode 0", "ExpecMode 1", "ExpecMode 2"
   :widths: 10, 20, 20, 20

   "256", "0.38", "0.11", "0.06"
   "1024", "2.44", "0.21", "0.19"
   "4096", "26.00", "1.67", "1.37"
   "16384", "424.5", "26.40", "20.84"

.. figure:: ../../../figs/fulldiag_expecmode_bench.png
   :name: fig_fulldiag_expecmode_bench_en
   :alt: Observable-evaluation time versus matrix dimension
   :width: 600px

   Observable-evaluation time versus matrix dimension :math:`N`
   (log-log).

At :math:`N=16384`, ``ExpecMode 1`` is about 16 times faster than
``ExpecMode 0`` (scaling with the process count), and ``ExpecMode 2``
gains another factor of about 1.3 from the trace kernels.

Multiple nodes
~~~~~~~~~~~~~~

Total run time for the 8-site Hubbard chain (:math:`N=4900`) with one-
and two-body Green-function output for all states, ``Solver 3`` (CPU),
on 2 nodes :math:`\times` 4 processes (AMD EPYC 7763, InfiniBand,
Intel MPI):

.. csv-table::
   :header: "ExpecMode", "Run time (2 nodes, 8 processes)"
   :widths: 20, 30

   "0 (conventional)", "11 min 52 s"
   "1 (state-task parallel)", "52 s"
   "2 (trace kernels)", "49 s"

The ``ExpecMode 0`` run time is dominated by the per-state inter-node
collectives (eigenvector gather and reductions :math:`\times N`
states); ``ExpecMode 1``/``2`` remove exactly this cost (about 14x
here). On multiple nodes we strongly recommend ``ExpecMode 1`` or
higher.

Guidelines
~~~~~~~~~~

* For small systems that fit in one process, ``Solver 0`` (LAPACK) is
  the simplest choice.
* From a few thousand dimensions upward with several processes,
  ``Solver 3`` (ELPA) wins on both diagonalization speed and memory
  (:math:`O(N^2/P)` via distributed Hamiltonian generation). If GPUs
  are available, ``NGPU`` accelerates it further. Note that the
  ``ExpecMode 1``/``2`` observable-evaluation speedups are also
  available with distributed ``Solver 1`` runs.
* For observables, use ``ExpecMode 2`` for the supported models and
  ``ExpecMode 1`` otherwise (the results agree with ``ExpecMode 0``).

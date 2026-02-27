MPI Parallelization
===================

This section describes the MPI parallelization strategy used in :math:`{\mathcal H}\Phi`.

Overview
--------

:math:`{\mathcal H}\Phi` uses MPI to parallelize calculations across multiple processes.
The Hilbert space basis states are distributed among MPI processes,
and Hamiltonian matrix-vector multiplication :math:`\hat{H}|\psi\rangle` is performed in parallel.

Site Classification
-------------------

When using :math:`N_{\rm proc}` MPI processes for an :math:`N`-site system,
sites are classified into two categories:

**Local sites** (:math:`N_{\rm local}` sites)
  Sites whose quantum states are stored within each process.
  For these sites, basis state manipulations can be performed locally.

**Inter-process sites** (:math:`N_{\rm inter}` sites)
  Sites whose quantum states span multiple processes.
  Operations on these sites require MPI communication.

The relationship is:

.. math::
   N_{\rm local} = N - \log_2(N_{\rm proc})

.. math::
   N_{\rm inter} = \log_2(N_{\rm proc})

For example, with 4 MPI processes (:math:`N_{\rm proc}=4`) on an 8-site system,
we have :math:`N_{\rm local}=6` local sites and :math:`N_{\rm inter}=2` inter-process sites.

MPI Communication Patterns
--------------------------

Operations in the Hamiltonian are classified by the MPI communication they require:

**Local operations**
  Both sites involved in the operation are local sites.
  No MPI communication is needed.

**MPIsingle operations**
  One site is local, one site is inter-process.
  Requires MPI communication with one partner process.
  The communication partner (origin) is determined by:

  .. math::
     \text{origin} = \text{myrank} \oplus T_{\rm pow}[\text{site}]

  where :math:`\oplus` denotes XOR operation and :math:`T_{\rm pow}` is the power of 2 for each site.

**MPIdouble operations**
  Both sites are inter-process sites.
  Requires MPI communication, but the received data can be processed uniformly.
  The communication partner is:

  .. math::
     \text{origin} = \text{myrank} \oplus (T_{\rm pow}[\text{site}_1] + T_{\rm pow}[\text{site}_2])

Batched MPI Communication
-------------------------

To reduce MPI communication overhead, :math:`{\mathcal H}\Phi` groups multiple operations
that share the same communication partner (origin) into batches:

**Before optimization:**
  Each transfer/interaction term with inter-process sites calls ``MPI_Sendrecv`` individually.
  If :math:`N` terms share the same origin, :math:`N` separate communications occur.

**After optimization:**
  Terms with the same origin are grouped together.
  A single ``MPI_Sendrecv`` retrieves the required data,
  then all :math:`N` terms are processed locally.

This optimization is particularly effective for:

* 2D/3D lattice models where many hopping terms share the same communication partner
* Systems with many transfer terms between neighboring sites

The batched communication is implemented for:

* SpinlessFermion models (MPIsingle)
* Hubbard models (MPIsingle, MPIdouble, InterAll)
* Spin models (Exchange interactions)

Fermion Sign
------------

For fermionic systems (Hubbard, SpinlessFermion), the fermion sign arising from
operator anticommutation must be carefully handled during MPI communication.

The sign is computed using the ``SgnBit`` function, which counts the number of
occupied states between the creation and annihilation operators.
This sign is pre-computed during initialization and stored with each transfer group.

Process Number Requirements
---------------------------

The number of MPI processes must satisfy specific constraints:

**Hubbard/Kondo models:**
  Process number = :math:`4^n` (due to 4 states per site: empty, up, down, double)

**Spin-1/2 models:**
  Process number = :math:`2^n` (due to 2 states per site: up, down)

**General spin models:**
  Process number = product of :math:`(2S_i+1)` for inter-process sites

See :ref:`Subsec:CreatingExpert` for detailed instructions on setting process numbers.

.. highlight:: none

.. _Subsec:TransSym:

TransSym file
-------------

This file restricts the calculation to one symmetry sector. It defines a
group :math:`G` of site permutations and a one-dimensional character
:math:`\chi(g)` (a complex number of unit modulus) for every element
:math:`g`. HPhi then works in the basis of symmetry-adapted states

.. math::

   |r;\chi\rangle \propto \sum_{g\in G}\chi(g)\,T_g|r\rangle ,

where :math:`T_g` moves the content of each site :math:`i` to the site
:math:`g(i)` and :math:`|r\rangle` runs over the representative
configurations. Every state of the sector satisfies
:math:`T_g|\psi\rangle=\chi(g)^{*}|\psi\rangle`, and the dimension of the
Hilbert space is reduced by roughly the order of the group. The sector
dimension is printed in the log as ``Symmetry basis: raw_dim=... sector_dim=...``.

The group is not restricted to translations. Reflections, rotations, and
their products with translations can be given in the same way, as long as
the Hamiltonian is invariant under every operation and the representation
is one-dimensional (abelian groups, or a one-dimensional irreducible
representation of a larger group). For a two-dimensional irreducible
representation, restrict the group to an abelian subgroup.

Standard mode writes this file as ``qptransidx.def`` when ``MomentumIndex``
is given (see :doc:`the MomentumIndex parameter <../standardmode_en/Parameters_for_conserved_quantities_en>`).
For a chain of length :math:`L` it contains the :math:`L` translations with
:math:`\chi(T^{g})=\exp(-2\pi i m g/L)`, which the ``MomentumIndex`` section
labels as the momentum :math:`k=2\pi m/L`.

An example of the file format (translations of a 6-site ring with
``MomentumIndex = 1``) is as follows.

::

    # MomentumIndex 1
    =============================================
    NQPTrans          6
    =============================================
    ======== TrIdx_TrWeight_and_TrIdx_i_xi ======
    =============================================
    0  1.000000000000000  0.000000000000000
    1  0.500000000000000 -0.866025403784439
    2 -0.500000000000000 -0.866025403784439
    3 -1.000000000000000  0.000000000000000
    4 -0.500000000000000  0.866025403784439
    5  0.500000000000000  0.866025403784439
    0 0 0 1
    0 1 1 1
    ...
    1 0 1 1
    1 1 2 1
    ...
    5 5 4 1

The following example is the reflection :math:`i\to 5-i` of the same ring,
with the character :math:`-1` (odd parity).

::

    =============================================
    NQPTrans          2
    =============================================
    ======== TrIdx_TrWeight_and_TrIdx_i_xi ======
    =============================================
    0  1.0
    1 -1.0
    0 0 0 1
    0 1 1 1
    0 2 2 1
    0 3 3 1
    0 4 4 1
    0 5 5 1
    1 0 5 1
    1 1 4 1
    1 2 3 1
    1 3 2 1
    1 4 1 1
    1 5 0 1

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Next [int01] lines: [int02] [double01] [double02]

*  Next [int01] :math:`\times` ``Nsite`` lines: [int02] [int03] [int04] [int05]

Lines that start with ``#`` and empty lines are skipped wherever they
appear and do not count in the line numbers above.

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the number of symmetry operations.
   Specify ``NQPTrans`` (case-insensitive).

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** The number of symmetry operations, i.e., the order of
   the group.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the index of a symmetry operation
   (:math:`0<=` [int02] :math:`<` [int01]).

*  [double01], [double02]

   **Type :** Double ([double02] can be omitted)

   **Description :** The real and imaginary parts of the character
   :math:`\chi(g)` of the operation [int02]. When [double02] is omitted,
   the imaginary part is zero.

*  [int03]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** A site index (:math:`0<=` [int03] :math:`<` ``Nsite``).

*  [int04]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** The site to which the operation [int02] moves the
   site [int03] (:math:`0<=` [int04] :math:`<` ``Nsite``).

*  [int05]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** Reserved for anti-periodic boundary conditions.
   In this version it must be ``1``.

Metadata
~~~~~~~~

A comment line of the form ``# MomentumIndex`` [int06] records the
``MomentumIndex`` that generated the file. [int06] must be a non-negative
integer representable by the C ``int`` type (at most ``INT_MAX``).
Standard mode writes it as the first line, and HPhi reports it in
the log as ``TransSym metadata: MomentumIndex=``\ [int06]. The value is not
used by the calculation in this version. Other comment lines are ignored.

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  Every operation must be a bijection of the sites, and every pair of
   [int02] and [int03] must be given exactly once.

*  The operations must form a group: they must contain the identity,
   and the composition of any two operations must be one of them.
   The characters must have unit modulus and be multiplicative,
   :math:`\chi(gh)=\chi(g)\chi(h)`, with :math:`\chi(e)=1` for the identity.

*  The Hamiltonian must be invariant under every operation. Each term is
   mapped by the operation and compared with the original terms; a
   mismatch terminates the program with
   ``TransSym Hamiltonian invariance failed``.

*  For ``SpinlessFermion`` and ``Hubbard``, the sign of the fermion
   permutation is taken into account automatically. For example, a
   reflection that exchanges occupied orbitals contributes a factor
   :math:`-1`, so the dimensions of the even and odd sectors differ from
   the counting for spins.

*  In this version the symmetry sector is available for ``Spin`` with
   :math:`S=1/2` and fixed ``2Sz`` (``Exchange`` and ``Ising`` terms),
   ``SpinlessFermion`` with fixed ``Ncond`` (``Trans`` and ``CoulombInter``
   terms), and ``Hubbard`` with fixed ``Nup`` and ``Ndown`` (``Trans`` and
   ``CoulombIntra`` terms), with the ``Lanczos`` and ``CG`` methods.
   Correlation functions, spectrum calculations, restart, and the input and
   output of Hamiltonians and eigenvectors are not supported together with
   this file. Unsupported combinations terminate the program with an error
   message that names the unsupported option.

.. raw:: latex

   \newpage

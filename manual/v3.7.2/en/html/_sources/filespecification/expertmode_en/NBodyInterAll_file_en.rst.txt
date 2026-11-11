.. highlight:: none

.. _Subsec:nbodyinterall:

NBodyInterAll file
------------------

This file determines generic N-body interaction terms added to the
Hamiltonian. A line with :math:`N` factors represents

.. math::

   {\mathcal H} += V \prod_{p=1}^{N}
   c_{i_p\sigma'_p}^{\dagger} c_{j_p\sigma_p}.

For Spin/SpinGC models, each factor represents a local spin matrix
element on one site instead of a fermion bilinear, and the input must
satisfy :math:`i_p=j_p`. The factors are written in the order used by
the operator product. An example of the file format is as follows.

::

    ==========================
    NNBodyInterAll 3
    ==========================
    ========NBodyInterAll=====
    ==========================
    1 0 0 1 0 -1.0 0.0
    1 1 0 0 0 -1.0 0.0
    2 0 0 0 0 1 1 1 1 0.5 0.0

.. _file_format_nbodyinterall:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:
   [int02] ([int03] [int04] [int05] [int06]) ... [double01] [double02].
   The group in parentheses is repeated [int02] times.

.. _parameters_nbodyinterall:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of generic N-body
   interaction terms. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of generic
   N-body interaction terms.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** The number of factors :math:`N` in this term.
   It must be a positive integer.

*  [int03], [int05] of each factor

   **Type :** Int (a blank parameter is not allowed)

   **Description :** Site indices
   (:math:`0<=` site index :math:`<` ``Nsite``). [int03] and [int05]
   correspond to :math:`i_p` and :math:`j_p`, respectively.

*  [int04], [int06] of each factor

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** Spin or local-state indices.
   | For Hubbard, tJ, Kondo, and their GC/NConserved variants:
   | 0: Up-spin,
   | 1: Down-spin.
   | For Spin/SpinGC, use the same local-state index convention as
     OneBodyG and TwoBodyG.
   | For SpinlessFermion/SpinlessFermionGC, only 0 is allowed
     (any non-zero index causes an error).

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for the real part of :math:`V`.

*  [double02]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for the imaginary part of :math:`V`.

.. _use_rules_nbodyinterall:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when [int01] is different from the total
   number of generic N-body interaction terms defined in this file.

*  A program is terminated when a site index or spin index is outside
   the allowed range.

*  For Spin/SpinGC, every factor must satisfy :math:`i_p=j_p`.

*  For Kondo/KondoGC/KondoNConserved, a factor involving a local-spin
   site must satisfy :math:`i_p=j_p`. General-spin Kondo local spins are
   not supported.

*  Diagonal NBodyInterAll terms must have a zero imaginary part.

*  Since the Hamiltonian must be Hermitian, off-diagonal NBodyInterAll
   terms must be written as adjacent Hermite-conjugate pairs. For
   fermion models, the conjugate line is written by reversing the factor
   order, swapping creation and annihilation indices in each factor, and
   using the complex-conjugate coefficient. For Spin/SpinGC, the strict
   Hermite-conjugate form, obtained by reversing the factor order,
   swapping the local-state indices, and using the complex-conjugate
   coefficient, is accepted. If all factors act on different sites, the
   spin operators commute, so the equivalent order-preserved form is also
   accepted.

   For example, in a fermion model, the following two adjacent lines form
   a Hermite-conjugate pair:

   ::

       2 0 0 1 0 2 1 3 1 0.25  0.10
       2 3 1 2 1 1 0 0 0 0.25 -0.10

   These lines correspond to

   .. math::

      \begin{aligned}
      &(0.25+0.10\mathrm{i})
      c_{0,0}^{\dagger} c_{1,0}
      c_{2,1}^{\dagger} c_{3,1},\\
      &(0.25-0.10\mathrm{i})
      c_{3,1}^{\dagger} c_{2,1}
      c_{1,0}^{\dagger} c_{0,0}.
      \end{aligned}

   In Spin/SpinGC, the strict Hermite-conjugate form can be written in the
   same way:

   ::

       2 0 1 0 0 2 0 2 1 0.25  0.10
       2 2 1 2 0 0 0 0 1 0.25 -0.10

   With :math:`X_i^{ab}=|a\rangle_i\langle b|`, these lines correspond to
   the following operators. For spin-1/2, the same pair can also be written
   with :math:`S^\pm`.

   .. math::

      \begin{aligned}
      &(0.25+0.10\mathrm{i}) X_0^{1,0} X_2^{0,1}
        =(0.25+0.10\mathrm{i}) S_0^+ S_2^-,\\
      &(0.25-0.10\mathrm{i}) X_2^{1,0} X_0^{0,1}
        =(0.25-0.10\mathrm{i}) S_2^+ S_0^-.
      \end{aligned}

   Since the two factors act on different sites, the second line may also
   keep the same site order as the first line:

   ::

       2 0 1 0 0 2 0 2 1 0.25  0.10
       2 0 0 0 1 2 1 2 0 0.25 -0.10

   This corresponds to the same Hermite-conjugate pair:

   .. math::

      \begin{aligned}
      &(0.25-0.10\mathrm{i}) X_0^{0,1} X_2^{1,0}
        =(0.25-0.10\mathrm{i}) S_0^- S_2^+
        =(0.25-0.10\mathrm{i}) S_2^+ S_0^-.
      \end{aligned}

*  NBodyInterAll is supported for SpinGC, Spin, HubbardGC, Hubbard,
   HubbardNConserved, SpinlessFermionGC, SpinlessFermion, tJGC, tJ,
   tJNConserved, KondoGC, Kondo, and KondoNConserved.

*  NBodyInterAll is not supported in TimeEvolution.

*  NBodyInterAll does not support CalcSpec for HubbardNConserved,
   tJNConserved, or KondoNConserved.

*  NBodyInterAll is not supported in FullDiag for
   SpinlessFermion/SpinlessFermionGC.

.. raw:: latex

   \newpage

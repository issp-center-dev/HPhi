.. highlight:: none

.. _Subsec:nbodyg:

NBodyG file
-----------

This file determines the target components of generic N-body Green's
functions. A line with :math:`N` factors represents

.. math::

   \left\langle \prod_{p=1}^{N}
   c_{i_p\sigma'_p}^{\dagger} c_{j_p\sigma_p} \right\rangle .

For Spin/SpinGC models, each factor represents a local spin matrix
element on one site instead of a fermion bilinear, and the input must
satisfy :math:`i_p=j_p`. An example of the file format is as follows.

::

    ==========================
    NNBodyG 3
    ==========================
    ============NBodyG========
    ==========================
    1 0 0 1 0
    2 0 0 0 0 1 1 1 1
    3 0 0 0 0 1 1 1 1 2 0 2 0

.. _file_format_nbodyg:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:
   [int02] ([int03] [int04] [int05] [int06]) ...
   The group in parentheses is repeated [int02] times.

.. _parameters_nbodyg:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of generic N-body
   Green's functions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of generic N-body
   Green's functions.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** The number of factors :math:`N` in this component.
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

.. _use_rules_nbodyg:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when [int01] is different from the total
   number of generic N-body Green's functions defined in this file.

*  A program is terminated when a site index or spin index is outside
   the allowed range.

*  For Spin/SpinGC, every factor must satisfy :math:`i_p=j_p`.

*  For Kondo/KondoGC/KondoNConserved, a factor involving a local-spin
   site must satisfy :math:`i_p=j_p`. General-spin Kondo local spins are
   not supported.

*  NBodyG is supported for SpinGC, Spin, HubbardGC, Hubbard,
   HubbardNConserved, SpinlessFermionGC, SpinlessFermion, tJGC, tJ,
   tJNConserved, KondoGC, Kondo, and KondoNConserved.

*  NBodyG is an equal-time expectation-value output and is not used as
   an excitation-operator input for CalcSpec. Use SingleExcitation or
   PairExcitation for dynamical correlation functions in CalcSpec.

*  In HubbardNConserved, tJNConserved, and KondoNConserved, using NBodyG
   together with CalcSpec is rejected.

*  The output values are written to the files described in
   :ref:`Subsec:nbodygdat`.

.. raw:: latex

   \newpage

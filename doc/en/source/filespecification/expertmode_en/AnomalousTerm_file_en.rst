.. highlight:: none

.. _Subsec:anomalousterm:

AnomalousTerm file
------------------

This file specifies anomalous pair Hamiltonian terms for HubbardGC
calculations. A line with type ``0`` represents
:math:`c_{i\sigma_1}c_{j\sigma_2}`, and a line with type ``1`` represents
:math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}`. An example of the
file format is as follows.

::

    ========================
    NAnomalousTerm 2
    ========================
    ========AnomalousTerm===
    ========================
    0 0 0 0 1 0.3000000000 0.0000000000
    1 0 1 0 0 0.3000000000 0.0000000000

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:
   [int02] [int03] [int04] [int05] [int06] [double01] [double02].

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of anomalous pair
   Hamiltonian terms. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of anomalous pair
   Hamiltonian terms.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** The anomalous pair type:
   | 0: :math:`c_{i\sigma_1}c_{j\sigma_2}`
   | 1: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}`.

*  [int03], [int05]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int03], [int05] :math:`<` ``Nsite``).

*  [int04], [int06]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | 0: Up-spin,
   | 1: Down-spin.

*  [double01], [double02]

   **Type :** Double (a blank parameter is not allowed)

   | **Description :** The coupling constant of the anomalous pair term.
   | [double01] and [double02] show the real and imaginary parts,
     respectively.

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  This file is currently supported only for HubbardGC calculations.

*  MPI execution is not supported for ``AnomalousTerm``.

*  ``AnomalousTerm`` is not supported for time evolution calculations or
   ``CalcSpec`` calculations.

*  The same fermion operator cannot be used twice in one anomalous pair.

*  Terms must appear as adjacent Hermitian conjugate pairs. For a term
   ``type i spin_i j spin_j V``, the next term must have the opposite type,
   the reversed operator order, and the complex conjugate coefficient:
   ``1-type j spin_j i spin_i conj(V)``.

*  A program is terminated when [int01] is different from the total number
   of anomalous pair Hamiltonian terms defined in this file.

*  A program is terminated when [int02]-[int06] are outside the range of the
   defined values.

.. raw:: latex

   \newpage

.. highlight:: none

.. _Subsec:pairexcitation:

PairExcitation file
-------------------

The operators to generate the pair excited state
:math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}(c_{i\sigma_1}^{\dagger}c_{j\sigma_2})`
are defined. The type of pair excitation operators
(:math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}` or
:math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`) must be same in the input
file. In the :math:`S_z` conserved system, :math:`\sigma_1` must be
equal to :math:`\sigma_2`. An example of the file format is as follows.

::

    ===============================
    NPair         24
    ===============================
    ======== Pair Excitation ======
    ===============================
        0     0     0     0    0    1.0    0.0
        0     1     0     1    0    1.0    0.0
        1     0     1     0    0    1.0    0.0
       (continue...)
       11     0    11     0    0    1.0    0.0
       11     1    11     1    0    1.0    0.0

.. _file_format_16:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:
   [int02]  [int03]  [int04]  [int05]  [int06]  [double01]  [double02].

.. _parameters_16:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of the pair
   excitation operators. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of pair
   excitation operators.

*  [int02], [int04]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int04] :math:`<` ``Nsite``).

*  [int03], [int05]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | For Hubbard or Kondo system, the index can be selected from
   | 0: Up-spin,
   | 1: Down-spin.
   | For Spin system, the index can be selected from
   | :math:`0, 1, \cdots, 2S+1` (corresponding to
     -:math:`S-0.5, -S+0.5, \cdots S+0.5`).

*  [int06]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a type of pair excitation
     operators:
   | 0: :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}`
   | 1: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`.

*  [double01], [double02]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** [double01] gives the real part
   of
   :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger} ( c_{i\sigma_1}^{\dagger}c_{j\sigma_2})`,
   while [double02] gives the imaginary part of
   :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger} ( c_{i\sigma_1}^{\dagger}c_{j\sigma_2})`.

.. _use_rules_16:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the pair excitation
   operators are double counted.

*  A program is terminated when the conditions :math:`i=j` and
   :math:`k=l` are not satisfied for the spin model.

*  A program is terminated when [int01] is different
   from the total number of two-body Green’s functions defined in this
   file.

*  A program is terminated when
   [int02]-[int06] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
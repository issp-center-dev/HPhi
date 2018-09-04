.. highlight:: none

.. _Subsec:singleexcitation:

SingleExcitation file
---------------------

The operators to generate the single excited state
:math:`c_{i\sigma_1}(c_{i\sigma_1}^{\dagger})` are defined. An example
of the file format is as follows.

::

    ===============================
    NSingle         24
    ===============================
    ======== Single Excitation ======
    ===============================
        0     0     0    1.0    0.0
        0     1     0    1.0    0.0
        1     0     0    1.0    0.0
       (continue...)
       11     0    0    1.0    0.0
       11     1    0    1.0    0.0

.. _file_format_15:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02]  [int03]  [int04]  [double01]  [double02].

.. _parameters_15:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of single excitation
   operators. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of total number
   of single excitation operators.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02] :math:`<` ``Nsite``).

*  [int03]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | For Hubbard or Kondo system, the index can be selected from
   | 0: Up-spin,
   | 1: Down-spin.
   | For Spin system, the index can be selected from
   | :math:`0, 1, \cdots, 2S+1` (corresponding to
     -:math:`S-0.5, -S+0.5, \cdots S+0.5`).

*  [int04]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a type of single excitation
     operators:
   | 0: :math:`c_{i\sigma_1}`
   | 1: :math:`c_{i\sigma_1}^{\dagger}`.

*  [double01], [double02]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** [double01] gives the real part
   of :math:`c_{i\sigma_1} (c_{i\sigma_1} ^{\dagger})`, while
   [double02] gives the imaginary part of
   :math:`c_{i\sigma_1} (c_{i\sigma_1} ^{\dagger})`.

.. _use_rules_15:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the single excitation
   operators are double counted.

*  A program is terminated when the conditions :math:`i=j` and
   :math:`k=l` are not satisfied for spin model.

*  A program is terminated when [int01] is different
   from the total number of two-body Green’s functions defined in this
   file.

*  A program is terminated, when
   [int02]-[int04] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
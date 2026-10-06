.. highlight:: none

.. _Subsec:twobodyg:

TwoBodyG file
-------------

This file determines the target components of the two-body Green’s
function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`.
For the spin system, the conditions :math:`i=j` and :math:`k=l` must be
satisfied. An example of the file format is as follows.

::

    =============================================
    NCisAjsCktAltDC        576
    =============================================
    ======== Green functions for Sq AND Nq ======
    =============================================
        0     0     0     0     0     0     0     0
        0     0     0     0     0     1     0     1
        0     0     0     0     1     0     1     0
        0     0     0     0     1     1     1     1
        0     0     0     0     2     0     2     0
        0     0     0     0     2     1     2     1
        0     0     0     0     3     0     3     0
        0     0     0     0     3     1     3     1
        0     0     0     0     4     0     4     0
        0     0     0     0     4     1     4     1
        0     0     0     0     5     0     5     0
        0     0     0     0     5     1     5     1
        0     0     0     0     6     0     6     0
        0     0     0     0     6     1     6     1
        0     0     0     0     7     0     7     0
        0     0     0     0     7     1     7     1
        0     0     0     0     8     0     8     0
        0     0     0     0     8     1     8     1
        0     0     0     0     9     0     9     0
        0     0     0     0     9     1     9     1
        0     0     0     0    10     0    10     0
        0     0     0     0    10     1    10     1
        0     0     0     0    11     0    11     0
        0     0     0     0    11     1    11     1
        0     1     0     1     0     0     0     0
        ...

.. _file_format_14:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:
   [int02]  [int03]  [int04]  [int05]  [int06]  [int07]  [int08]  [int09].

.. _parameters_14:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of two-body Green’s
   functions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of two-body
   Green’s functions.

*  [int02],
   [int04],[int06],
   [int08]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int04], [int06], [int08] :math:`<` ``Nsite``).

*  [int03],
   [int05],[int07],
   [int09]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | For Hubbard or Kondo system, the index can be selected from
   | 0: Up-spin,
   | 1: Down-spin.
   | For Spin system, the index can be selected from
   | :math:`0, 1, \cdots, 2S+1` (corresponding to
     -:math:`S-0.5, -S+0.5, \cdots S+0.5`).

.. _use_rules_14:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the two-body Green’s
   functions are double counted.

*  A program is terminated when the conditions :math:`i=j` and
   :math:`k=l` are not satisfied for the spin model.

*  A program is terminated when [int01] is different
   from the total number of two-body Green’s functions defined in this
   file.

*  A program is terminated, when
   [int02]-[int09] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
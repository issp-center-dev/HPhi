.. highlight:: none

.. _Subsec:onebodyg:

OneBodyG file
-------------

This file determines the target components of the one-body Green’s
function :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`. An
example of the file format is as follows.

::

    ===============================
    NCisAjs         24
    ===============================
    ======== Green functions ======
    ===============================
        0     0     0     0
        0     1     0     1
        1     0     1     0
        1     1     1     1
        2     0     2     0
        2     1     2     1
        3     0     3     0
        3     1     3     1
        4     0     4     0
        4     1     4     1
        5     0     5     0
        5     1     5     1
        6     0     6     0
        6     1     6     1
        7     0     7     0
        7     1     7     1
        8     0     8     0
        8     1     8     1
        9     0     9     0
        9     1     9     1
       10     0    10     0
       10     1    10     1
       11     0    11     0
       11     1    11     1

.. _file_format_13:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02]  [int03]  [int04]  [int05].

.. _parameters_13:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of one-body Green’s
   functions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of one-body
   Green’s functions.

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

.. _use_rules_13:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the one-body Green’s
   functions are double counted.

*  A program is terminated when [int01] is different
   from the total number of one-body Green’s functions defined in this
   file.

*  A program is terminated when
   [int02]-[int05] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
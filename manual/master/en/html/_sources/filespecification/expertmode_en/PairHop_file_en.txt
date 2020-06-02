.. highlight:: none

PairHop file
------------

This file determines the values of PairHop couplings
:math:`J_{ij}^{\rm Pair}` (for :math:`S=1/2` system),

.. math::

   {\mathcal H}+=\sum_{i,j}J_{ij}^{\rm Pair} (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}+h.c.).

An example of the file format is as follows.

::

    ====================== 
    NPairhop 6 
    ====================== 
    ========Pairhop ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

.. _file_format_9:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [int03] [double01].

.. _parameters_9:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of PairHop
   couplings. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of PairHop
   couplings.

*  [int02], [int03]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int03] :math:`<` ``Nsite``).

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm Pair}`.

.. _use_rules_9:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when [int01] is different
   from the total number of the PairHop couplings defined in this file.

*  A program is terminated when either [int02] or
   [int03] is outside the range of the defined
   values.

.. raw:: latex

   \newpage
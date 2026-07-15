.. highlight:: none

.. _Subsec:pairlift:

PairLift file
-------------

This file determines the values of PairLift couplings
:math:`J_{ij}^{\rm PairLift}` (for :math:`S=1/2` system),

.. math::

   {\mathcal H}+=\sum_{i,j}J_{ij}^{\rm PairLift} (c_ {i \uparrow}^{\dagger}c_{i\downarrow}c_{j \uparrow}^{\dagger}c_{j \downarrow}+c_ {i \downarrow}^{\dagger}c_{i\uparrow}c_{j \downarrow}^{\dagger}c_{j \uparrow}).

An example of the file format is as follows.

::

    ====================== 
    NPairLift 6  
    ====================== 
    ========NPairLift ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

.. _file_format_12:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [int03] [double01] .

.. _parameters_12:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of PairLift
   couplings. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of PairLift
   couplings.

*  [int02], [int03]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int03] :math:`<` ``Nsite``).

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm PairLift}`.

.. _use_rules_12:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the PairLift couplings
   are double counted.

*  A program is terminated when [int01] is different
   from the total number of PairLift couplings defined in this file.

*  A program is terminated when either [int02] or
   [int03] is outside the range of the defined
   values.

.. raw:: latex

   \newpage
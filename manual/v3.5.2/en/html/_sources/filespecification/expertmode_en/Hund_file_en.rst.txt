.. highlight:: none

Hund file
---------

This file determines the values of Hund couplings
:math:`J_{ij}^{\rm Hund}` (for :math:`S=1/2` system),

.. math::

   \mathcal H+=-\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}).

An example of the file format is as follows.

::

    ====================== 
    NHund 6  
    ====================== 
    ========Hund ====== 
    ====================== 
       0     1 -0.250000
       1     2 -0.250000
       2     3 -0.250000
       3     4 -0.250000
       4     5 -0.250000
       5     0 -0.250000

.. _file_format_8:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [int03] [double01].

.. _parameters_8:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of Hund couplings.
   You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of Hund
   couplings.

*  [int02], [int03]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int03] :math:`<` ``Nsite``).

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm Hund}`.

.. _use_rules_8:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the Hund couplings are
   double counted.

*  A program is terminated when [int01] is different
   from the total number of Hund couplings defined in this file.

*  A program is terminated when either [int02] or
   [int03] is outside the range of the defined
   values.

.. raw:: latex

   \newpage
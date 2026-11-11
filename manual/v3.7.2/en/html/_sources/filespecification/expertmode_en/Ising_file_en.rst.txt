.. highlight:: none

Ising file
----------

This file determines the values of Ising interactions :math:`J_{ij}^{z}`
(for :math:`S=1/2` system). For the fermion electronic system, the Ising
terms are given as

.. math::

   {\mathcal H}+=\sum_{i,j}J_{ij}^{z} (n_{i\uparrow}-n_{i\downarrow})(n_{j\uparrow}-n_{j\downarrow} ).

For the spin system, they are given as

.. math::

   {\mathcal H}+=\sum_{i,j}J_{ij}^{z} S_ {i}^{z}S_{j}^z.

An example of the file format is as follows.

::

    ====================== 
    NIsing 6  
    ====================== 
    ========Ising ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

.. _file_format_11:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [int03] [double01].

.. _parameters_11:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of Ising
   interactions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of Ising
   interactions.

*  [int02], [int03]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int03] :math:`<` ``Nsite``).

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for :math:`J_{ij}^{\rm z}`.

.. _use_rules_11:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the Ising interactions
   are double counted.

*  A program is terminated when [int01] is different
   from the total number of Ising interactions defined in this file.

*  A program is terminated when either [int02] or
   [int03] is outside the range of the defined
   values.

.. raw:: latex

   \newpage
.. highlight:: none

CoulombIntra file
-----------------

This file determines the values of the on-site interactions :math:`U_i`
(for :math:`S=1/2` system),

.. math::

   {\mathcal H}+=\sum_{i}U_i n_ {i \uparrow}n_{i \downarrow}.

An example of the file format is as follows.

::

    ====================== 
    NCoulombIntra 6  
    ====================== 
    ========i_0LocSpn_1IteElc ====== 
    ====================== 
       0  4.000000
       1  4.000000
       2  4.000000
       3  4.000000
       4  4.000000
       5  4.000000

.. _file_format_6:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [double01].

.. _parameters_6:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of on-site
   interactions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of on-site
   interactions.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02] :math:`<` ``Nsite``).

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for :math:`U_i`.

.. _use_rules_6:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of on-site interactions
   are double counted.

*  A program is terminated when [int01] is different
   from the total number of on-site interactions defined in this file.

*  A program is terminated when [int02] is outside
   the range of the defined values.

.. raw:: latex

   \newpage
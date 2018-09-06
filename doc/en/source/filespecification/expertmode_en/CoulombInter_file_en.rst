.. highlight:: none

CoulombInter file
~~~~~~~~~~~~~~~~~

This file determines the values of off-site interactions :math:`V_{ij}`
(for :math:`S=1/2` system),

.. math:: H+=\sum_{i,j}V_{ij} n_ {i}n_{j}.

 An example of the file format is as follows.

::

    ====================== 
    NCoulombInter 6  
    ====================== 
    ========CoulombInter ====== 
    ====================== 
       0     1  1.0000
       1     2  1.0000
       2     3  1.0000
       3     4  1.0000
       4     5  1.0000
       5     0  1.0000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: [string01] [int01]

-  Lines 3-5: Header

-  Lines 6-: [int02] [int03] [double01].

Parameters
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of off-site
   interactions. You can freely give a name to the keyword.

-  :math:`[`\ int01\ :math:`]`

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of off-site
   interactions.

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int03\ :math:`]`

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<= [`\ int02\ :math:`], [`\ int03\ :math:`]<\verb|Nsite|`).

-  :math:`[`\ double01\ :math:`]`

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for :math:`V_{ij}`.

Use rules
^^^^^^^^^

-  Headers cannot be omitted.

-  A program is terminated when the components of off-site interactions
   are double counted.

-  A program is terminated when :math:`[`\ int01\ :math:`]` is different
   from the total number of off-site interactions defined in this file.

-  A program is terminated when either :math:`[`\ int02\ :math:`]` or
   :math:`[`\ int03\ :math:`]` is outside the range of the defined
   values.

.. raw:: latex

   \newpage

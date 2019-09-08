.. highlight:: none

.. _Subsec:Trans:

Trans file
----------

This file determines the values of the transfer integrals :math:`t_{ij\sigma_1\sigma2}`,

.. math::

   {\mathcal H} +=-\sum_{ij\sigma_1\sigma2} t_{ij\sigma_1\sigma2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}.

An example of the file format is as follows.

::

    ======================== 
    NTransfer      24  
    ======================== 
    ========i_j_s_tijs====== 
    ======================== 
        0     0     2     0   1.000000  0.000000
        2     0     0     0   1.000000  0.000000
        0     1     2     1   1.000000  0.000000
        2     1     0     1   1.000000  0.000000
        2     0     4     0   1.000000  0.000000
        4     0     2     0   1.000000  0.000000
        2     1     4     1   1.000000  0.000000
        4     1     2     1   1.000000  0.000000
        4     0     6     0   1.000000  0.000000
        6     0     4     0   1.000000  0.000000
        4     1     6     1   1.000000  0.000000
        6     1     4     1   1.000000  0.000000
        6     0     8     0   1.000000  0.000000
        8     0     6     0   1.000000  0.000000
    ...

.. _file_format_4:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02]  [int03]  [int04]  [int05]  [double01]  [double02].

.. _parameters_4:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of transfer
   integrals. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of transfer
   integrals.

*  [int02], [int04]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int04] :math:`<` ``Nsite`` ).

*  [int03], [int05]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | 0: Up-spin
   | 1: Down-spin.

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for a real part of
   :math:`t_{ij\sigma_1\sigma_2}`.

*  [double02]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for an imaginary part of
   :math:`t_{ij\sigma_1\sigma_2}`.

.. _use_rules_4:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  Since the Hamiltonian must be Hermitian, the relation
   :math:`t_{ij\sigma_1\sigma_2}=t_{ji\sigma_2\sigma_1}^{\dagger}` must
   be satisfied. A program is terminated when this relation is broken.

*  A program is terminated when the components of the on-site
   interactions are double counted.

*  A program is terminated when [int01] is different
   from the total number of transfer integrals defined in this file.

*  A program is terminated when
   [int02]-[int05] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
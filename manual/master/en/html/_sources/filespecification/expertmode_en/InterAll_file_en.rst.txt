.. highlight:: none

.. _Subsec:interall:

InterAll file
-------------

This file determines the values of generalized two-body interactions
integrals :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`,

.. math::

   {\mathcal H}+=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
   I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}.

For spin, the conditions :math:`i=j` and :math:`k=l` must be satisfied. An example of the file format is as follows.

::

    ====================== 
    NInterAll      36  
    ====================== 
    ========zInterAll===== 
    ====================== 
    0    0    0    1    1    1    1    0   0.50  0.0
    0    1    0    0    1    0    1    1   0.50  0.0
    0    0    0    0    1    0    1    0   0.25  0.0
    0    0    0    0    1    1    1    1  -0.25  0.0
    0    1    0    1    1    0    1    0  -0.25  0.0
    0    1    0    1    1    1    1    1   0.25  0.0
    2    0    2    1    3    1    3    0   0.50  0.0
    2    1    2    0    3    0    3    1   0.50  0.0
    2    0    2    0    3    0    3    0   0.25  0.0
    2    0    2    0    3    1    3    1  -0.25  0.0
    2    1    2    1    3    0    3    0  -0.25  0.0
    2    1    2    1    3    1    3    1   0.25  0.0
    4    0    4    1    5    1    5    0   0.50  0.0
    4    1    4    0    5    0    5    1   0.50  0.0
    4    0    4    0    5    0    5    0   0.25  0.0
    4    0    4    0    5    1    5    1  -0.25  0.0
    4    1    4    1    5    0    5    0  -0.25  0.0
    4    1    4    1    5    1    5    1   0.25  0.0
    ...

.. _file_format_5:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:
   [int02] [int03] [int04] [int05] [int06] [int07] [int08] [int09] [double01] [double02].

.. _parameters_5:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of generalized
   two-body interactions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of generalized
   two-body interactions.

*  [int02], [int04],
   [int06], [int08]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int04], [int06], [int08] :math:`<` ``Nsite``).

*  [int03], [int05],
   [int07], [int09]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | 0: Up-spin,
   | 1: Down-spin.

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for a real part of
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`.

*  [double02]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for an imaginary part of
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`.

.. _use_rules_5:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  Since the Hamiltonian must be Hermitian, the relation
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}=I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}^{\dagger}`
   must be satisfied. A program is terminated when this relation is
   broken. It is noted that the term of the Hermitian conjugate for
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}`
   should be inputted as
   :math:`I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}`
   :math:`c_{l\sigma_4}^{\dagger}c_{k\sigma_3}c_{j\sigma_2}^{\dagger}c_{i\sigma_1}`.

*  A program is terminated when the conditions :math:`i=j` and
   :math:`k=l` are not satisfied for the spin model.

*  A program is terminated when the components of the on-site
   interactions are double counted.

*  A program is terminated when [int01] is different
   from the total number of generalized two-body interactions defined in
   this file.

*  A program is terminated when
   [int02]-[int09] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
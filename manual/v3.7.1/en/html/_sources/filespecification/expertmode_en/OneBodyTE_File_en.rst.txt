.. highlight:: none

.. _Subsec:1TE:

OneBodyTE File
--------------

This file determines the values of the transfer integrals :math:`t_{ij\sigma_1\sigma_2} (t)` at each time :math:`t`,

.. math::

   {\mathcal H}(t) +=-\sum_{ij\sigma_1\sigma_2} t_{ij\sigma_1\sigma_2}(t)c_{i\sigma_1}^{\dagger}c_{j\sigma_2}.

An example of the file format is as follows.

::

    ==================================
    AllTimeStep         100
    ==================================
    ===== OneBody Time Evolution =====
    ==================================
        0.0   4
        0  0  1  0      1.0    0.0
        1  0  0  0      1.0    0.0
        0  1  1  1      1.0    0.0
        1  1  0  1      1.0    0.0
        0.2   4
       (continue...)

.. _file_format_18:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

From line 6, time :math:`t` and total number of transfer integrals
:math:`N(t)` are first defined and next the transfer integrals
:math:`t_{ij\sigma_1\sigma_2} (t)` are defined.

*  Line m: [double01]  [int02]

*  Lines (m+1) - (m+1+[int02]):
   [int03]  [int04]  [int05]  [int06]  [double02]  [double03]

.. _parameters_18:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of transfer
   integrals. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total time steps defined in
   this file.

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** Time :math:`t`.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of transfer
   integrals at time :math:`t`.

*  [int03], [int05]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int03], [int05] :math:`<` ``Nsite``).

*  [int04], [int06]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | 0: Up-spin
   | 1: Down-spin.

*  [double02], [double03]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for real and imaginary part of
   :math:`t_{ij\sigma_1\sigma_2}(t)` at time :math:`t` is defined
   [double02] and [double03],
   respectively.

.. _use_rules_17:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when ``LanczosStep`` defined in ``ModPara``
   is greater than [int02].

*  A program is terminated when
   [int03]-[int06] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
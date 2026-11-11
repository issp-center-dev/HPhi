.. highlight:: none

.. _Subsec:2TE:

TwoBodyTE File
---------------

This file determines the values of the two-body interactions :math:`I_{i\sigma_1j \sigma_2 k \sigma_3 l \sigma_4}(t)c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}` at each time :math:`t`,

.. math::

   {\mathcal H}(t) +=\sum_{ijkl}\sum_{\sigma_1\sigma_2\sigma_3\sigma_4} I_{i\sigma_1j \sigma_2 k \sigma_3 l \sigma_4}(t) c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}.

An example of the file format is as follows.

::

    ==================================
    AllTimeStep         100
    ==================================
    ===== TwoBody Time Evolution =====
    ==================================
        0.0   3
        0  1  0  1  0  0  0  0       1.0    0.0
        1  0  1  0  0  0  0  0       1.0    0.0
        1  0  1  0  2  0  2  0       1.0    0.0 
        0.2   3
       (continue...)

.. _file_format_19:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

From line 6, time :math:`t` and total number of transfer integrals
:math:`N(t)` are first defined and next the two-body interactions
:math:`I_{i\sigma_1j \sigma_2 k \sigma_3 l \sigma_4}(t)` are defined.

*  Line m: [double01]  [int02]

*  Lines (m+1) - (m+1+[int02]):
   [int03]  [int04]  [int05]  [int06]  [int07]  [int08]  [int09]  [int10]  [double02]  [double03]

.. _parameters_19:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of generalized
   two-body interactions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total time steps defined in
   this file.

*  [double01]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** Time :math:`t`.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of two-body
   interactions at time :math:`t`.

*  [int03], [int05],
   [int07], [int09]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int03],[int05],[int07],[int09] :math:`<` ``Nsite``).

*  [int04], [int06],
   [int08], [int10]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | 0: Up-spin
   | 1: Down-spin.

*  [double02], [double03]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** A value for real and imaginary part of
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}(t)` at time :math:`t`
   is defined [double02] and
   [double03], respectively.

.. _use_rules_18:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when ``LanczosStep`` defined in ``ModPara``
   is greater than [int02].

*  A program is terminated when
   [int03]-[int10] are outside
   the range of the defined values.

.. raw:: latex

   \newpage
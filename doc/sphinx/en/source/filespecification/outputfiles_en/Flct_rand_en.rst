.. highlight:: none

.. _Subsec:flctrand:

Flct_rand.dat
-------------

| (For the TPQ method) This file is outputted to show the calculation
  results of the fluctuation of the particle number, doublon, and
  :math:`S_z` for the TPQ method. In the restart calculation, the values
  are added to the previous file. An example of the file format is as
  follows.

::

     # inv_temp, N, N^2, D, D^2, Sz, Sz^2, step_i
    0.0826564 12.00 144.00 0.00 0.00 0.0009345626081113 0.2500 1
    0.1639935 12.00 144.00 0.00 0.00 0.0023147006319775 0.2500 2
    0.2440168 12.00 144.00 0.00 0.00 0.0037424057659867 0.2500 3
    ...
    135.97669 12.00 144.00 0.00 0.00 -0.0000000000167368 0.2500 1998
    136.04474 12.00 144.00 0.00 0.00 -0.0000000000165344 0.2500 1999

.. _file_name_9:

File name
~~~~~~~~~

*  Flct_rand??.dat

?? indicates the number of runs in the calculation of the TPQ method.

.. _file_format_32:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [double03]
   [double04] [double05]
   [double06] [double07]
   [int01].

.. _parameters_32:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** Inverse temperature :math:`1/{k_{\rm B}T}`.

*  [double02]

   **Type :** Double

   **Description :** A total particle number
   :math:`\sum_{i} \langle \hat{n}_i \rangle`.

*  [double03]

   **Type :** Double

   **Description :** The expected value of the square of the particle
   number :math:`\langle (\sum_{i} \hat{n}_i)^2 \rangle`.

*  [double04]

   **Type :** Double

   **Description :** The expected value of doublon
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   (:math:`N_s` is the total number of sites).

*  [double05]

   **Type :** Double

   **Description :** The expected value of the square of doublon
   :math:`\frac{1}{N_s}\langle ( \sum_{i} n_{i\uparrow} n_{i\downarrow})^2\rangle`
   (:math:`N_s` is the total number of sites).

*  [double06]

   **Type :** Double

   **Description :** The expected value of :math:`S_z`
   :math:`\frac{1}{N_s} \sum_{i}\langle \hat{S}_i^z\rangle` (:math:`N_s`
   is the total number of sites).

*  [double07]

   **Type :** Double

   **Description :** The expected value of the square of :math:`S_z`
   :math:`\frac{1}{N_s} \langle (\sum_{i} \hat{S}_i^z)^2\rangle`
   (:math:`N_s` is the total number of sites).

*  [int01]

   **Type :** Int

   **Description :** The number of operations of
   :math:`(l-\hat{\mathcal H}/N_{s})` for an initial wave function, where
   :math:`l` is ``LargeValue`` defined in a ModPara file and
   :math:`N_{s}` is the total number of sites.

.. raw:: latex

   \newpage
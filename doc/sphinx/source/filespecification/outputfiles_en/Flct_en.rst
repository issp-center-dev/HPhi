.. highlight:: none

.. _Subsec:flct:

Flct.dat
--------

| (For the real time evolution method) This file is outputted to show
  the calculation results of the fluctuation of the particle number,
  doublon, and :math:`S_z`. In the restart calculation, the values are
  added to the previous file. An example of the file format is as
  follows.

::

     # time, N, N^2, D, D^2, Sz, Sz^2, step_i
    0.0000000000000000 7.9999999999999991 63.9999999999999929 ...
    0.0100000000000000 8.0000000000000604 64.0000000000004832 ...
    0.0200000000000000 8.0000000000000018 64.0000000000000142 ...
    0.0300000000000000 8.0000000000001013 64.0000000000008100 ...
    0.0400000000000000 7.9999999999999183 63.9999999999993463 ...
    0.0500000000000000 7.9999999999999520 63.9999999999996163 ...
    0.0600000000000000 7.9999999999999627 63.9999999999997016 ...
    0.0700000000000000 8.0000000000000835 64.0000000000006679 ...
    0.0800000000000000 8.0000000000000924 64.0000000000007390 ...
    0.0900000000000000 7.9999999999999600 63.9999999999996803 ...
    0.1000000000000000 7.9999999999999067 63.9999999999992539 ...
    ...

.. _file_name_13:

File name
~~~~~~~~~

*  Flct.dat

.. _file_format_36:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [double03]
   [double04] [double05]
   [double06] [double07]
   [int01].

.. _parameters_36:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** Time :math:`t`

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

   **Description :** Time step.

.. raw:: latex

   \newpage
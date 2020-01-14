.. highlight:: none

.. _Subsec:ssrand:

SS_rand.dat
-----------

| (For the TPQ method) This file is outputted to show the calculation
  results for the TPQ method. In the restart calculation, the values are
  added to the previous file. An example of the file format is as
  follows.

::

    # inv_tmp, energy, phys_var, phys_doublon, phys_num, step_i
    0.017471  5.526334 45.390269 1.464589 6.000000 1
    0.034863  5.266718 42.655559 1.434679 6.000000 2
    ...
    31.999572  -4.814170 23.176231 0.590568 6.000000 1997
    32.015596  -4.814170 23.176231 0.590568 6.000000 1998
    32.031620  -4.814170 23.176231 0.590568 6.000000 1999

.. _file_name_8:

File name
~~~~~~~~~

*  SS_rand??.dat

?? indicates the number of runs in the calculation of the TPQ method.

.. _file_format_31:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [double03]
   [double04] [double05]
   [int01].

.. _parameters_31:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** Inverse temperature :math:`1/{k_{\rm B}T} (k_{\rm B} = 1)`.

*  [double02]

   **Type :** Double

   **Description :** The expected value of the energy
   :math:`\langle \mathcal H \rangle`.

*  [double03]

   **Type :** Double

   **Description :** The expected value of the square of the Hamiltonian
   :math:`\langle \mathcal H^2 \rangle`.

*  [double03]

   **Type :** Double

   **Description :** The expected value of the doublon,
   :math:`\sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`.

*  [double05]

   **Type :** Double

   **Description :** The total number of particles
   :math:`\langle {\hat n} \rangle`.

*  [int01]

   **Type :** Int

   **Description :** The number of operations of
   :math:`(l-\hat{\mathcal H}/N_{s})` for an initial wave function, where
   :math:`l` is ``LargeValue`` defined in a ModPara file and
   :math:`N_{s}` is the total number of sites.

.. raw:: latex

   \newpage
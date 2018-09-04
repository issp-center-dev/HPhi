.. highlight:: none

.. _Subsec:normrand:

Norm_rand.dat
-------------

| (For the TPQ method) This file is outputted to show the calculation
  process information for the TPQ method. In the restart calculation,
  the values are added to the previous file. An example of the file
  format is as follows.

::

     # inv_temp, global_norm, global_1st_norm, step_i 
    0.017471 19.046586 11.288975 1
    0.034863 19.089752 11.288975 2
    ...
    31.999572 20.802362 11.288975 1997
    32.015596 20.802362 11.288975 1998
    32.031620 20.802362 11.288975 1999

.. _file_name_7:

File name
~~~~~~~~~

*  Norm_rand??.dat

?? indicates the number of runs in the calculation of the TPQ method.

.. _file_format_30:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [double03]
   [int01].

.. _parameters_30:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** Inverse temperature :math:`1/{k_{\rm B}T}`.

*  [double02]

   **Type :** Double

   **Description :** The norm of a wave function before normalization
   given by :math:`\langle \tilde{\psi}_{k} |\tilde{\psi}_{k}\rangle`,
   where
   :math:`|\tilde{\psi}_{k}\rangle \equiv(l-\hat{\mathcal H}/N_{s})|\psi_{k=1}\rangle`.

*  [double03]

   **Type :** Double

   **Description :** The norm of an initial wave function before
   normalization given by
   :math:`\langle \tilde{\psi}_{0} |\tilde{\psi}_{0}\rangle`, where
   :math:`|\tilde{\psi}_{0}\rangle` is an initial random vector.

*  [int01]

   **Type :** Int

   **Description :** The number of operations of
   :math:`(l-\hat{\mathcal H}/N_{s})` for an initial wave function, where
   :math:`l` is ``LargeValue`` defined in a ModPara file and
   :math:`N_{s}` is the total number of sites.

.. raw:: latex

   \newpage
.. highlight:: none

.. _Subsec:ssTE:

SS.dat
------

| (For the real time evolution method) This file is outputted to show
  the calculation results. In the restart calculation, the values are
  added to the previous file. An example of the file format is as
  follows.

::

     # time, energy, phys_var, phys_doublon, phys_num, step_i
    0.0000000000000000  -6.0412438187293001 38.8635272290786489 ...
    0.0100000000000000  -5.9310482979751606 37.9593669819686923 ...
    0.0200000000000000  -5.8287182777288828 37.1390062697724801 ...
    0.0300000000000000  -5.7384311863736031 36.4282319427381651 ...
    0.0400000000000000  -5.6636677030535481 35.8466140292489897 ...
    0.0500000000000000  -5.6070659264425551 35.4081795274108799 ...
    0.0600000000000000  -5.5703150294725914 35.1222606981659666 ...
    0.0700000000000000  -5.5540895348193438 34.9942503380419154 ...
    0.0800000000000000  -5.5580244678717312 35.0260574979670665 ...
    0.0900000000000000  -5.5807312981978212 35.2161499389042660 ...
    0.1000000000000000  -5.6198544688797947 35.5591788356216298 ...
    ...

.. _file_name_12:

File name
~~~~~~~~~

*  SS.dat

.. _file_format_35:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [double03]
   [double04] [double05]
   [int01].

.. _parameters_35:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** Time :math:`t`.

*  [double02]

   **Type :** Double

   **Description :** The expected value of the energy
   :math:`\langle H \rangle`.

*  [double03]

   **Type :** Double

   **Description :** The expected value of the square of the Hamiltonian
   :math:`\langle H^2 \rangle`.

*  [double04]

   **Type :** Double

   **Description :** The expected value of the doublon,
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   (:math:`N_{s}` is the total number of sites).

*  [double05]

   **Type :** Double

   **Description :** The total number of particles
   :math:`\langle {\hat n} \rangle`.

*  [int01]

   **Type :** Int

   **Description :** Time step.

.. raw:: latex

   \newpage
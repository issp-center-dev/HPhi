.. highlight:: none

.. _subsec:energy.dat:

energy.dat
----------

| (For the Lanczos method) The values of the energy, doublon, and
  :math:`\langle S_z \rangle` calculated by using the eigenvector
  obtained by the Lanczos or CG method are outputted. An example of the
  file format is as follows.
| For ``method="Lanczos"``

::

    Energy  -7.1043675920 
    Doublon  0.4164356536 
    Sz  0.0000000000 

For ``method="CG"``

::

    State 0
      Energy  -7.1043675920 
      Doublon  0.4164356536 
      Sz  0.0000000000 

    State 1
    :

.. _file_name_4:

File name
~~~~~~~~~

*  ##_energy.dat

## indicates a header defined by [string02] in a ModPara file.

.. _file_format_27:

File format
~~~~~~~~~~~

*  Line 1: Energy [double01]

*  Line 2: Doublon [double02]

*  Line 3: Sz [double02].

.. _parameters_27:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** The value of the energy calculated by the
   eigenvetor obtained by the Lanczos or CG method.

*  [double02]

   **Type :** Double

   **Description :** The value of the doublon calculated by the
   eigenvetor obtained by the Lanczos or CG method,
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   (:math:`N_s` is the total number of sites).

*  [double03]

   **Type :** Double

   **Description :** The value of :math:`S_z` calculated by the
   eigenvetor obtained by the Lanczos or CG method.

.. raw:: latex

   \newpage



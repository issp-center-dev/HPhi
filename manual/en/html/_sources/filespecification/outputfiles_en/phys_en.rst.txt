.. highlight:: none

.. _Subsec:phys:

phys.dat
--------

| (For the FullDiag method) This file is outputted to show the physical
  values calculated by the FullDiag method. The data are outputted in
  descending order for energies. An example of the file format is as
  follows.

::

     <H>         <N>        <Sz>       <S2>       <D> 
      -4.814170   0.000000   0.000000  -0.000000   0.590568
      -3.796850   0.000000   0.000000   1.333333   0.423804
     ...
     14.489622   0.000000   0.000000   0.000000   2.550240
     14.852520   0.000000   0.000000   0.000000   2.329157

.. _file_name_14:

File name
~~~~~~~~~

*  Canonical ensemble: ##_phys_Nup_$$Ndown%%.dat

*  Grand canonical ensemble: ##_phys.dat.

##, $$, and %% indicate [string02], Nup, and Ndown defined in a ModPara
file, respectively.

.. _file_format_38:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [double03]
   [double04] [double05].

.. _parameters_38:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** The energy :math:`\langle \mathcal H\rangle`.

*  [double02]

   **Type :** Double

   **Description :** The total number of particles
   :math:`\langle \hat{n}\rangle`.

*  [double03]

   **Type :** Double

   **Description :** The expected value of :math:`S_z`,
   :math:`\langle S_z\rangle`.

*  [double04]

   **Type :** Double

   **Description :** The expected value of :math:`{\boldsymbol S}^2` ,
   :math:`\langle {\boldsymbol S}^2 \rangle`.

*  [double05]

   **Type :** Double

   **Description :** The expected value of doublon,
   :math:`\frac{1}{N_s} \sum_{i}\langle n_{i\uparrow}n_{i\downarrow}\rangle`
   (:math:`N_{s}` is the total number of sites).

.. raw:: latex

   \newpage
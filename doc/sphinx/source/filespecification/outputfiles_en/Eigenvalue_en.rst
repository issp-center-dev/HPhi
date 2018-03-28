.. highlight:: none

.. _Subsec:eigenvalue:

Eigenvalue.dat
--------------

| (For the FullDiag method) This file is outputted to show the energies
  calculated by the FullDiag method. An example of the file format is as
  follows.

::

     0 -4.8141698096 
     1 -3.7968502453 
     2 -3.2462822372 
     ...
     397 13.9898305290 
     398 14.4896221034 
     399 14.8525199079 

.. _file_format_37:

File format
~~~~~~~~~~~

*  [int01] [double01]

.. _parameters_37:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The index of eigenvalues. The index 0 is for the
   energy of the ground state and the indexes are labeled in descending
   order for energies.

*  [double01]

   **Type :** Double

   **Description :** The expected value of energy
   :math:`\langle \mathcal H \rangle`.

.. raw:: latex

   \newpage
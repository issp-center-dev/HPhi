.. highlight:: none

WarningOnTransfer.dat
---------------------

This file shows the double counted components of transfer integrals. An
example of the file format is as follows.

::

    double conuntings in transfers: i=0 j=2 spni 0 spnj 0  
    double conuntings in transfers: i=2 j=0 spni 0 spnj 0  
    double conuntings in transfers: i=0 j=2 spni 1 spnj 1  
    double conuntings in transfers: i=2 j=0 spni 1 spnj 1  

.. _file_format_26:

File format
~~~~~~~~~~~

*  Double countings in transfers: i=[int01]
   j=[int02] spni [int03] spnj 
   [int04]

.. _parameters_26:

Parameters
~~~~~~~~~~

*  [int01], [int02]

   **Type :** Int

   **Description :** The integer of the site number where the transfer
   integrals are double counted.

*  [int03], [int04]

   **Type :** Int

   | **Description :** The integer of the spin index of a transfer
     integral:
   | 0: Up-spin
   | 1: Down-spin.

.. raw:: latex

   \newpage

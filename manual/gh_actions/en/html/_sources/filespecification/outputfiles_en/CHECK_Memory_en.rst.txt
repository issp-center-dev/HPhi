.. highlight:: none

CHECK_Memory.dat
----------------

This file shows the size of the memory used in the calculation. An
example of the file format is as follows.

::

    MAX DIMENSION idim_max=400 
    REQUIRED MEMORY  max_mem=0.000019 GB 

.. _file_format_25:

File format
~~~~~~~~~~~

*  MAX DIMENSION idim_max=[int01]

*  REQUIRED MEMORY max_mem =[double01] GB

.. _parameters_25:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** An integer to show the total numbers of the Hilbert
   space under a calculation.

*  [double01]

   **Type :** Double

   **Description :** The size of memory to store the Hilbert space in a
   calculation (GB unit).

.. raw:: latex

   \newpage
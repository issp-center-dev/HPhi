.. highlight:: none

.. _Subsec:normTE:

Norm.dat
--------

| (For the real time evolution method) This file is outputted to show
  the calculation process information for the real time evolution
  method. In the restart calculation, the values are added to the
  previous file. An example of the file format is as follows.

::

     # time, global_norm, global_1st_norm, step_i
    0.0000000000000000 0.9999999999999994 0.0000000000000000 0
    0.0100000000000000 1.0000233421898765 0.0000000000000000 1
    0.0200000000000000 1.0000211100654208 0.0000000000000000 2
    0.0300000000000000 1.0000182214014706 0.0000000000000000 3
    0.0400000000000000 1.0000148460317946 0.0000000000000000 4
    0.0500000000000000 1.0000111372562455 0.0000000000000000 5
    0.0600000000000000 1.0000072252313270 0.0000000000000000 6
    0.0700000000000000 1.0000032174168609 0.0000000000000000 7
    0.0800000000000000 0.9999992048382456 0.0000000000000000 8
    0.0900000000000000 0.9999952720176869 0.0000000000000000 9
    0.1000000000000000 0.9999915078951970 0.0000000000000000 10

.. _file_name_11:

File name
~~~~~~~~~

*  Norm.dat

.. _file_format_34:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Lines 2-: [double01]
   [double02] [int01].

.. _parameters_34:

Parameters
~~~~~~~~~~

*  [double01]

   **Type :** Double

   **Description :** Time :math:`t`

*  [double02]

   **Type :** Double

   **Description :** The norm of a wave function at :math:`t`.

*  [int01]

   **Type :** Int

   **Description :** Time step.

.. raw:: latex

   \newpage
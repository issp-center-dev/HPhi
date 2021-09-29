.. highlight:: none

Time_TPQ_Step.dat
-----------------

| (For the TPQ method) This file is outputted to show the time for
  starting the calculation of the TPQ method at each seed and step. In
  the restart calculation, the values are added to the previous file. An
  example is as follows.

::

    set 0 step 1:TPQ begins: Wed Jul 13 07:59:20 2016
    set 0 step 2:TPQ begins: Wed Jul 13 07:59:20 2016
    set 0 step 3:TPQ begins: Wed Jul 13 07:59:20 2016
    ...
    set 4 step 1997:TPQ begins: Wed Jul 13 07:59:32 2016
    set 4 step 1998:TPQ begins: Wed Jul 13 07:59:32 2016
    set 4 step 1999:TPQ begins: Wed Jul 13 07:59:32 2016

.. _file_name_6:

File name
~~~~~~~~~

*  ##_TPQ_Step.dat

## indicates a header defined by [string02] in a ModPara file.

.. _file_format_29:

File format
~~~~~~~~~~~

*  Set [int01] stp [int02]: TPQ
   begins: [string01]

.. _parameters_29:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The seed number in the calculation of the TPQ
   method.

*  [int02]

   **Type :** Int

   **Description :** The step number in the calculation of the TPQ
   method.

*  [string01]

   **Type :** String

   **Description :** The time for starting the calculation of the TPQ
   method at each seed and step.

.. raw:: latex

   \newpage
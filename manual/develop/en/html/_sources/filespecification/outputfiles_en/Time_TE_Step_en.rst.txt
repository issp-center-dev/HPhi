.. highlight:: none

Time_TE_Step.dat
----------------

| (For the real time evolution method) This file is outputted to show
  the time for starting the calculation of the real time evolution
  method at each step. In the restart calculation, the values are added
  to the previous file. An example is as follows.

::

    step 1:TE begins: Wed Jul 13 07:59:20 2017
    step 2:TE begins: Wed Jul 13 07:59:20 2017
    step 3:TE begins: Wed Jul 13 07:59:20 2017
    …
    step 1997:TE begins: Wed Jul 13 07:59:32 2017
    step 1998:TE begins: Wed Jul 13 07:59:32 2017
    step 1999:TE begins: Wed Jul 13 07:59:32 2017

.. _file_name_10:

File name
~~~~~~~~~

*  ##_TE_Step.dat

## indicates a header defined by [string02] in a ModPara file.

.. _file_format_33:

File format
~~~~~~~~~~~

*  stp [int01]: TE begins:
   [string01]

.. _parameters_33:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The step number in the calculation.

*  [string01]

   **Type :** String

   **Description :** The time for starting the calculation at each step.

.. raw:: latex

   \newpage
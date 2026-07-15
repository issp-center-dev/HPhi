.. highlight:: none

.. _Subsec:locspn:

LocSpin file
------------

| This file determines sites with localized spins. The file format is as
  follows.

::

    ================================ 
    NlocalSpin     6  
    ================================ 
    ========i_0LocSpn_1IteElc ====== 
    ================================ 
        0      1
        1      0
        2      1
        3      0
        4      1
        5      0
        6      1
        7      0
        8      1
        9      0
       10      1
       11      0

.. _file_format_3:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [int03].

.. _parameters_3:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of localized spins.
   You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of localized
   spins.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02] :math:`<` ``Nsite``).

*  [int03]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer for selecting an electron state
     whether the electron state is a localized spin or an itinerant
     electron state:
   |  0: Itinerant electron state
     :math:`n>`\0: Localized spin state with :math:`2S=n`.

.. _use_rules_3:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when [int01] is different
   from the total number of localized spins indicated by
   [int03].

*  A program is terminated, when [int02] is
   different from the total number of sites.

*  A program is terminated under the condition 
   [int02] :math:`<0` or 
   ``Nsite`` :math:`<=` [int02].

.. raw:: latex

   \newpage
.. highlight:: none

.. _Subsec:excitedvec:

eigenvec.dat
------------

When OutputExcitedVec=1 in a CalcMod file, the excited vector calculated by using the input vector and excited operators defined in parr.def is outputted.
An example of the file format is as follows.

::

   4096
   0.0009135367, -0.0004231751
   -0.0009712430, 0.0293999045
   -0.0294644178, 0.0210086435
   0.0214247977, -0.0107587147
   0.0272643022, -0.0404634256
   0.0034322694, -0.0184446640
   0.0019911098, 0.0004403706
   0.0114685735, -0.0114381935
   0.0092888239, 0.0088235535
        …

File name
~~~~~~~~~

*  ##_excitedvec\_&&\_rank\_$$.dat

## indicates [string02] in a ModPara file, && is the number of
eigenvalues, and $$ is a number of rank.


File format
~~~~~~~~~~~

*  Line 1:[int01] 

*  Lines2-:[double01] [double02]


Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** Dimenssions of the excited vector.

*  [double01], [double02]

   **Type :** Double

   **Description :** The value of the excited vector;
   [double01] and [double02]
   represent real and imaginary part of the excited vector, respectively.



   
.. raw:: latex

   \newpage

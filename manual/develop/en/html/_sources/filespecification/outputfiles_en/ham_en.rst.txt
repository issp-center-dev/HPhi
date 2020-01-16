.. highlight:: none

.. _Subsec:ham:

ham.dat
-------

| (For the FullDiag method) When ``OutputHam=1`` in the ``CalcMod``
  file, the Hamiltonian calculated by :math:`{\mathcal H}\Phi` is outputted by the
  MatrixMarket format. The recalculation by using this file can be down
  when ``InputHam=1`` in the ``CalcMod`` file. An example of the file
  format is as follows.

::

    %%%%MatrixMarket matrix coordinate complex hermitian
    28 28 56
    1 1 1.000000 0.000000
    2 1 0.500000 0.000000
    3 2 0.500000 0.000000
    4 3 0.500000 0.000000
    5 4 0.500000 0.000000
    6 5 0.500000 0.000000
    7 6 0.500000 0.000000
    7 7 1.000000 0.000000
        …

.. _file_name_15:

File name
~~~~~~~~~

*  ##_ham.dat

## indicates [string02] in a ModPara file.

.. _file_format_39:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2:
   [int01] [int02] [int03]

*  Lines3-:[int04] [int05] [double01] [double02]

.. _parameters_39:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** Number of rows of Hamiltonian.

*  [int02]

   **Type :** Int

   **Description :** Number of columns of Hamiltonian.

*  [int03]

   **Type :** Int

   **Description :** Number of nonzero elements of Hamiltonian.

*  [double01], [double02]

   **Type :** Double

   **Description :** The value of Hamiltonian;
   [double01] and [double02]
   represent real and imaginary part of the Hamiltonian, respectively.

.. raw:: latex

   \newpage
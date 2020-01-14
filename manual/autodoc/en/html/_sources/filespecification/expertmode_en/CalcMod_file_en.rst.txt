.. highlight:: none

.. _Subsec:calcmod:

CalcMod file
------------

This file determines the parameters for the calculation method, model, and output mode. The file format is as follows.

::

    CalcType   0
    CalcModel   2
    CalcEigenVec 0

.. _file_format_1:

File format
~~~~~~~~~~~

[string01] [int01]

.. _parameters_1:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String

   **Description :** Select a word from keywords.

*  [int01]

   **Type :** Int

   | **Description :** A parameter that is correlated with a keyword.

.. _use_rules_1:

Use rules
~~~~~~~~~

*  After setting the keywords at [string 01], a half-width blank is
   needed for setting a parameter.

*  Keywords can be set in random order.

*  If the keywords or filenames are incorrect, the program is
   terminated.

*  The keywords “CalcType" and “CalcModel" are essential.

*  When a head of line is \"#", the line is skipped.

 

Keywords and parameters
~~~~~~~~~~~~~~~~~~~~~~~

The parameters correlated with the keywords are as follows.

*  ``CalcType``

   **Type :** Int

   | **Description :** Select the method for calculation from the
     following list:
   | 0: Lanczos method
   | 1: Analysis of the physical properties by using TPQ
   | 2: Full diagonalization method.
   | 3: LOBCG for the ground state.

*  ``CalcModel``

   **Type :** Int

   | **Description :** Select the model from the following list:
   | 0: Fermion Hubbard model (canonical ensemble: conservation of
     particles or conservation of particles and the component of
     :math:`S_z`)
   | 1: Spin model (canonical ensemble: conservation of the component of
     :math:`S_z`)
   | 2: Kondo lattice model (canonical ensemble: conservation of
     particles, the component of :math:`S_z`)
   | 3: Fermion Hubbard model (grand canonical ensemble)
   | 4: Spin model (grand canonical ensemble)
   | 5: Kondo lattice model (grand canonical ensemble).

   For the fermion Hubbard model, you can select the model under the
   conservation of the particles by setting ``NCond`` in the ModPara
   file. When you want to select the model under the conservation of
   particles and the component of :math:`S_z`, set both ``NCond`` and
   ``2Sz`` in the ModPara file.

*  ``CalcEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the method to calculate the eigenvectors:
   | 0: Lanczos+CG methods (when the convergence of eigenvectors is not
     sufficient for using the Lanczos method, the CG method is applied
     to calculate eigenvectors).
   | 1: Lanczos method.

*  ``InitialVecType``

   **Type :** Int (default value: 0)

   | **Description :** Select the type of an initial vector:
   | 0: Complex type
   | 1: Real type.

*  ``OutputEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of outputting an eigenvector:
   | 0: Not output an eigenvector
   | 1: Output an eigenvector.

*  ``InputEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of inputting an eigenvector:
   | 0: Not input an eigenvector
   | 1: Input an eigenvector.

*  ``ReStart``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of inputting a restart vector:
   | 0: Not restart calculation
   | 1: Output a restart vector
   | 2: Input a restart vector and output a new restart vector
   | 3: Input a restart vector.

*  ``CalcSpec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of calculating dynamical Green’s functions:
   | 0: Not calculate dynamical Green’s functions
   | 1: (not restart) Input an initial vector and files for generating single excited or pair excited states
   | 2: Input components of triangular diagonal matrix
   | 3: Output both components of triangular diagonal matrix and a restart vector
   | 4: Input both components of triangular diagonal matrix and a restart vector
   | 5: Input and output both components of triangular diagonal matrix and a restart vector.

*  ``OutputHam``

   **Type :** Int (default value: 0)

   | **Description :** Full Diag)Select the mode of outputting Hamiltonian:
   | 0: not output Hamiltonian.
   | 1: output Hamiltonian.

*  ``InputHam``

   **Type :** Int (default value: 0)

   | **Description :** (Full Diag)Select the mode of inputting Hamiltonian:
   | 0: not input Hamiltonian.
   | 1: input Hamiltonian.

*  ``OutputExcitedVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the mode of outputting an excited vector:
   | 0: Not output an eigenvector
   | 1: Output an eigenvector.
   
*  ``Scalapack``

   **Type :** Int (default value: 0)

   | **Description :** (Full Diag)Select to use ScaLAPACK library for full diagonalization:
   | 0: not to use ScaLAPACK.
   | 1: use ScaLAPACK.


*  ``NGPU``

   **Type :** Int (default value: 2)

   | **Description :** (Full Diag)Select the number of GPU devices for full diagonalization:
   | :math:`{\mathcal H} \Phi` does not support to use GPU devices at multi-nodes. 

.. raw:: latex

   \newpage

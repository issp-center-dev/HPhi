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
   | 1: mTPQ method
   | 2: Full diagonalization method
   | 3: LOBCG for the ground state
   | 4: Time-evolution
   | 5: cTPQ method

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
   | 7: Spinless fermion model (canonical ensemble: conservation of particles)
   | 8: Spinless fermion model (grand canonical ensemble).
   | 9: :math:`t`-:math:`J` model (canonical ensemble: conservation of
     particles, or conservation of particles and the component of
     :math:`S_z`)
   | 10: :math:`t`-:math:`J` model (grand canonical ensemble)

   For the fermion Hubbard model, you can select the model under the
   conservation of the particles by setting ``NCond`` in the ModPara
   file. When you want to select the model under the conservation of
   particles and the component of :math:`S_z`, set both ``NCond`` and
   ``2Sz`` in the ModPara file.

   The :math:`t`-:math:`J` models (9, 10) follow the same ``NCond`` /
   ``2Sz`` selection as the fermion Hubbard model: for the canonical model
   (9), setting only ``NCond`` conserves the total number of electrons,
   while setting both ``NCond`` and ``2Sz`` also conserves :math:`S_z`.
   Doubly-occupied sites are excluded from the Hilbert space (the local
   dimension per site is 3: empty, up, or down), and the models are
   available only in the expert mode. Note that for MPI runs the number of
   processes must nevertheless be a power of four, the same as for the
   Fermion Hubbard model (**not** a power of three), because the internal
   representation keeps four states per site.

   For the spinless fermion model, only Trans (hopping) and CoulombInter
   (inter-site interaction) terms are valid. CoulombIntra, Hund, Exchange,
   and PairHop cannot be used since there are no spin degrees of freedom.

*  ``CalcEigenVec``

   **Type :** Int (default value: 0)

   | **Description :** Select the method to calculate the eigenvectors:
   | 0: Lanczos+CG methods (when the convergence of eigenvectors is not
     sufficient for using the Lanczos method, the CG method is applied
     to calculate eigenvectors).
   | 1: Lanczos method.

*  ``InitialVecType``

   **Type :** Int (default value: 0)

   | **Description :** Select the type of an initial vector (:math:`v0`):
   | -1: Real part (:math:`{\rm Re}[v0]]`) and imaginary part  (:math:`{\rm Re}[v0]]`) of the initial 
    vector are give as the normally distributed random numbers. Thus, the normalized initial vectors are uniformly distributed
    on the :math:`N_{\rm H}` dimensional super sphere (:math:`N_{\rm H}` is the dimension of the Hilbert space). 
   | 0: Complex type (:math:`{\rm Re}[v0]\in[-1:1]`, :math:`{\rm Im}[v0]\in[-1:1]` ).
   | 1: Real type (:math:`{\rm Re}[v0]\in[-1:1]`, :math:`{\rm Im}[v0]=0`).

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
   
*  ``OutputDataHead``

   **Type :** Int (default value: 0)

   | **Description :** Select whether to prefix TPQ/TE physical quantity output filenames (``SS``, ``Norm``, ``Flct``) with the header string defined by ``CDataFileHead`` in the ModPara file:
   | 0: Do not add a prefix (e.g., ``SS_rand0.dat``).
   | 1: Add the ``CDataFileHead`` prefix (e.g., ``zvo_SS_rand0.dat``).

*  ``OutputGreenFormat``

   **Type :** Int (default value: 0)

   | **Description :** Select the output format for Green function files:
   | 0: Existing split files.
   | 1: Aggregate indexed files for TPQ/cTPQ, real-time evolution, Full diagonalization, and LOBCG.
   | In aggregate mode, TPQ/cTPQ rows start with ``set`` and ``step``, real-time evolution rows start with ``step``, and Full diagonalization/LOBCG rows start with ``eigen``.
   | ``AnomalousG`` in LOBCG keeps the existing non-aggregate output because the existing output is not split by eigen index.

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

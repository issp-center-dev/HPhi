.. highlight:: none

.. _Subsec:modpara:

ModPara file
-----------------

| This file determines the parameters for calculation. The file format
  is as follows.

::

    --------------------
    Model_Parameters   0
    --------------------
    VMC_Cal_Parameters
    --------------------
    CDataFileHead  zvo
    CParaFileHead  zqp
    --------------------
    Nsite          16   
    Ncond          16    
    2Sz            0 
    Lanczos_max    1000 
    initial_iv     12   
    exct           1    
    LanczosEps     14   
    LanczosTarget  2    
    LargeValue     12   
    NumAve         5    
    ExpecInterval  20   

.. _file_format_2:

File format
~~~~~~~~~~~

*  Lines 1-4: Header

*  Line 6: [string01] [string02]

*  Lines 7-8: Header

*  Lines 9- : [string01] [int01].

.. _parameters_2:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String

   **Description :** Select a word from keywords.

*  [string02]

   **Type :** String (a blank parameter is not allowed)

   **Description :** Set a header for output files.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** A parameter that is correlated with a keyword.

.. _use_rules_2:

Use rules
~~~~~~~~~

*  From Line 9: After setting keywords at [string 01], a half-width
   blank is needed for setting a parameter

*  All the parameters are needed and the order for the parameters is
   fixed

 

.. _keywords_and_parameter_1:

Keywords and parameters
~~~~~~~~~~~~~~~~~~~~~~~

In the following, common parameters and parameters for each method are
shown.

 

Common parameters
~~~~~~~~~~~~~~~~~

*  ``CDataFileHead``

   **Type :** String (a blank parameter is not allowed)

   **Description :** A header for output files. For example, the output
   filename for one-body Green’s function becomes
   "**xxx_Lanczos_cisajs.dat**" (xxx are the characters set by
   ``CDataFileHead``).

*  ``Nsite``

   **Type :** Int (positive integer)

   **Description :** The number of sites.

*  ``Ncond``

   **Type :** Int (positive integer)

   **Description :** The number of conduction electrons (not used in
   grand canonical ensemble).

*  ``2Sz``

   **Type :** Int (positive integer)

   **Description :** The total value of :math:`2S_z` (not used in grand
   canonical ensemble). For conservation of :math:`S_z` in the case of
   ``CalcModel`` = 0 (fermion Hubbard model) or 2 (Kondo lattice model),
   set ``Ncond``.

*  ``initial_iv``

   **Type :** Int

   **Description :** An initial vector is specified with this parameter.

   *  Lanczos method

      *  For canonical ensemble and ``initial_iv`` :math:`\geq 0`

         The non-zero components of an initial vector are specified with
         this parameter.

      *  For grand canonical ensemble or ``initial_iv`` :math:`< 0`

         The seed of the random generator is given by this parameter and
         the random vector is used as the initial vector.

   *  TPQ method

      The seed of the random generator is given by this parameter and
      the random vector is used as the initial vector.

   See :ref:`Ch:algorithm` for details of setting an
   initial vector.

*  ``CalcHS``

   **Type :** Int (positive integer)

   **Description :** If CalcHS=1, an efficient algorithm for generating
   the restricted Hilbert space with the specified quantum number is
   used (Details of algorithm is shown in
   http://qlms.github.io/HPhi/develop/tips.pdf [in Japanese]). Default
   value is 1 and the efficient algorithm is used.

 

Lanczos method
~~~~~~~~~~~~~~

*  ``Lanczos_max``

   **Type :** Int (positive integer)

   **Description :** The maximum number of Lanczos steps in the
   calculation. When the convergence within the specified accuracy is
   satisfied, the calculation is completed before a step reaches
   ``Lanczos_max``. In the case of restart calculation, ``Lanczos_max``
   must be larger than that of the previous calculation.

*  ``exct``

   **Type :** Int (positive integer)

   **Description :** An integer for setting the number of eigenvectors
   obtained from the ground energy by the Lanczos method.

*  ``LanczosEps``

   **Type :** Int (positive integer)

   **Description :** An integer for judging the convergence of the
   Lanczos method. The convergence is determined by whether the
   condition is satisfied that the relative error between an eigenvalue
   and an eigenvalue at the Lanczos step of the one step before is less
   than :math:`10^{- \verb|LanczosEps|}`.

*  ``LanczosTarget``

   **Type :** Int (positive integer)

   **Description :** An integer giving the target of the eigenvalue for
   judging the convergence of the Lanczos method. For example, the
   target becomes a ground state when ``LanczosTarget`` is equal to one,
   and a first excited state when ``LanczosTarget`` is equal to two.

 

CG method
~~~~~~~~~

*  ``exct``

   **Type :** Int (positive integer)

   **Description :** The number of eigenvectors is specified.

*  ``Lanczos_max``

   **Type :** Int (positive integer)

   **Description :** The maximum number of iteration steps in the
   calculation. When the convergence within the specified accuracy is
   satisfied, the calculation is completed before a step reaches
   Lanczos_max. In the case of restart calculation, ``Lanczos_max``
   must be larger than that of the previous calculation.

*  ``LanczosEps``

   **Type :** Int (positive integer)

   **Description :** For ``method="CG"``, the calculation finishes when
   the 2-norm of the residual vector becomes smaller than
   :math:`10^{- \verb|LanczosEps|/2}`.

 

TPQ method
~~~~~~~~~~

*  ``Lanczos_max``

   **Type :** Int (positive integer)

   **Description :** The total number of TPQ steps is specified with
   this parameter. In the case of restart calculation, ``Lanczos_max``
   must be larger than that of the previous calculation.

*  ``LargeValue``

   **Type :** Double

   **Description :** An integer giving :math:`l` of 
   :math:`l=\hat{\mathcal H}/N_{s}` used in the TPQ method.

*  ``NumAve``

   **Type :** Int

   **Description :** An integer giving the number of independent runs
   for the TPQ method.

*  ``ExpecInterval``

   **Type :** Int

   | **Description :** An integer giving the interval steps of
     calculating the correlation functions in the TPQ method.
   | **Note:** A small interval increases the time cost of calculations.

 

Calculating dynamical Green’s functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``OmegaOrg``

   **Type :** Complex

   **Description :** The center value of the frequency. Specify the real
   and imaginary parts in that order separated by a space, and if there
   is no imaginary part, the real part of the frequency is only given.

*  ``OmegaIm``

   **Type :** Double

   **Description :** The imaginary part of the frequency. When
   ``OmegaOrg`` is defined in a ``modpara`` file, ``OmegaIm`` is added
   to the imaginary value of ``OmegaOrg``.

*  ``OmegaMin``

   **Type :** Complex

   **Description :** The lower limit of the frequency from ``OmegaOrg``.
   Specify the real and imaginary parts in that order separated by a
   space, and if there is no imaginary part, the real part of the
   frequency is only given.

*  ``OmegaMax``

   **Type :** Complex

   **Description :** The upper limit of the frequency from 
   ``OmegaOrg``. Specify the real and imaginary parts in that order
   separated by a space, and if there is no imaginary part, a real part
   of the frequency is only given.

*  ``NOmega``

   **Type :** Int

   **Description :** The integer for defining the step size of the 
   frequency :math:`\Delta \omega = (` ``OmegaMax`` - 
   ``OmegaMin`` :math:`)/N_{\omega}`. The frequency is given by
   :math:`z_n=` ``OmegaOrg``\ :math:`+`\ ``OmegaMin``\ :math:`+ \Delta \omega \times n`.

Real time evolution method
~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``Lanczos_max``

   **Type :** Int (positive integer)

   **Description :** The total number of real time evolution steps is
   specified with this parameter. In the case of restart calculation,
   ``Lanczos_max`` must be larger than that of the previous calculation.

*  ``ExpandCoef``

   **Type :** Int (positive integer)

   **Description :** An integer giving the expansion order :math:`n` for
   real time evolution method;

   .. math:: \exp\left(-i \hat{\cal H} \Delta t \right) = \sum_{i=0}^{N}\frac{1}{n!}\left(-i \hat{\cal H} \Delta t \right)^n .

*  ``ExpecInterval``

   **Type :** Int (positive integer)

   | **Description :** An integer giving the interval steps of
     calculating the correlation functions.
   | **Note:** A small interval increases the time cost of calculations.

*  ``OutputInterval``

   **Type :** Int (positive integer)

   | **Description :** An integer giving the interval steps of output
     the wave function.
   | The wave vector is output when ``OutputEigenVec=1`` in ``CalcMod``
     file.

.. raw:: latex

   \newpage

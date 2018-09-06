.. highlight:: none

Parameters for the numerical condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``2S``

   **Type :** Positive integer (``1`` as a default)

   **Description :** The :math:`2 S` at each site in the localized spin
   system is specified. (E.g. ``1`` for the :math:`1/2` system)

*  ``Restart``

   **Type :** String (choose from ``"None"``, ``"Restart_out"``,
   ``"Restart_in"``, ``"Restart"``. ``"None"`` as a default)

   **Description :** The condition of the restart is specified.
   ``"None"`` for omitting file IOs for the restart, ``"Restart_out"``
   for starting calculation from scratch and generating a restart-file
   after the calculation finishes, ``"Restart_in"`` for starting
   calculation with the restart-file generated in the previous run,
   ``"Restart"`` for ``"Restart_out"`` + ``"Restart_in"``.

*  ``anczos_max``

   **Type :** Positive integer (default value: ``2000``)

   **Description :** The upper limit of the Lanczos/LOBCG/BiCG step and
   the number of steps for TPQ/Time=evolution are specified with this
   parameter.

*  ``initial_iv``

   **Type :** Integer (default value: ``-1``)

   **Description :** An initial vector is specified with this parameter.

   *  Lanczos method

      *  For the canonical ensemble and ``initial_iv`` :math:`\geq 0`

         The non-zero components of an initial vector are specified with
         this parameter.

      *  For the grand canonical ensemble or ``initial_iv`` :math:`< 0`

         The seed of the random generator is given by this parameter and
         the random vector is used as the initial vector.

   *  TPQ method

      The seed of the random generator is given by this parameter and
      the random vector is used as the initial vector.

   See :ref:`Ch:Algorithm` for details of setting an
   initial vector.

*  ``exct``

   **Type :** Positive integer (default value: ``1``)

   | **Description :** The number of eigenvectors obtained from the
     ground energy by the Lanczos method are specified.
   | When exct=2, the eigenvector of the first-excited state is
     obtained. When ``method="CG"``, the number of states to be
     calculated is specified.

   **Note**: The condition ``nvec`` :math:`>=` ``exct`` must be
   satisfied.

*  ``LanczosEps``

   **Type :** Positive integer (default value: ``14``)

   **Description :** The convergence criterion for the Lanczos method is
   specified with this parameter. If the difference between the old and
   the new target eigenvalue falls below
   :math:`10^{- LanczosEps|}`, the Lanczos step will finish. For
   ``method="CG"``, we assume the calculation is converged when the
   2-norm of the residual vector becomes smaller than
   :math:`10^{-{\tt LanczosEps}/2}`.

*  ``LanczosTarget``

   **Type :** Positive integer (default value: ``2``)

   **Description :** The target eigenenergy for the convergence
   criterion is specified. If it is set to ``1``, the target eigenenergy
   becomes the ground state.

*  ``LargeValue``

   **Type :** Double (the default value is provided below)

   **Description :** (Only for TPQ) :math:`l` as :math:`l=\hat{\mathcal H}/N_{s}`
   is used in the TPQ calculation. Usually, the largest eigenvalue of
   the Hamiltonian is used as :math:`l`. Thus, the default value of
   :math:`l` is taken as the summation of the absolute values of each
   coefficient in the Hamiltonian divided by the number of sites.

*  ``NumAve``

   **Type :** Positive integer (default value: ``5``)

   **Description :** (Only for TPQ) The number of independent runs for
   the TPQ method is specified with this parameter.

*  ``ExpecInterval``

   **Type :** Positive integer (default value: ``20``)

   | **Description :** (Only for TPQ) The interval of calculating
     correlation functions in the TPQ iteration is specified.
   | **Note:** A small interval increases the time cost of calculations.

*  ``OutputMode``

   **Type :** Choose from ``"none"``, ``"correlation"``, and ``"full"``
   (``correlation`` as default)

   **Description :** Indices of correlation functions are specified with
   this keyword. ``"none"`` indicates correlation functions will not be
   calculated. When ``outputmode="correlation"``, the correlation
   function supported by the utility ``fourier`` is computed. For more
   details, see the document in ``doc/fourier/``. If ``"full"`` is
   selected, :math:`\langle c_{i \sigma}^{\dagger}c_{j \sigma'} \rangle`
   is computed at all :math:`i, j, \sigma, \sigma'`, and
   :math:`\langle c_{i_1 \sigma_1}^{\dagger}c_{i_2 \sigma_2} c_{i_3 \sigma_3}^{\dagger}c_{i_4 \sigma_4} \rangle`
   is computed at all
   :math:`i_1, i_2, i_3, i_4, \sigma_1, \sigma_2, \sigma_3, \sigma_4`.

   In a spin system, the indices are specified as those of the
   Bogoliubov representation (see :ref:`Sec:sec_bogoliubov_rep` ).

*  ``InitialVecType``

   **Type :** Character (choose from ``"C"``, ``"R"``. ``"C"`` as a
   default)

   **Description :** The type of the initial eigenvector is specified.
   ``C`` for the complex number, and ``R`` for the real number.

*  ``EigenVecIO``

   **Type :** String (choose from ``"None"``, ``"Out"``, ``"In"``.
   ``"None"`` as a default)

   **Description :** The I/O of the eigenvector is specified. ``"None"``
   for omitting the IO of the eigenvector, ``"Out"`` for writing the
   eigenvector to a file, ``"In"`` for reading the eigenvector from a
   file and using it in the subsequent calculation (such as the Green’s
   function).

.. raw:: latex

   \newpage

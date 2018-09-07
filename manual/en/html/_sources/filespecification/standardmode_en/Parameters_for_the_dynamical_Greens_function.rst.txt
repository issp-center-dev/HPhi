.. highlight:: none

Parameters for the dynamical Green’s function
---------------------------------------------

* ``CalcSpec``

   **Type :** String(choose from ``"None"``, ``"Normal"``,
   ``"NoIteration"``, ``"Restart_out"``, ``"Restart_in"``,
   ``"Restart"``. ``"None"`` as default.)

   **Description :** The condition for the calculation of the dynamical
   Green’s function is specified. ``"None"`` for omitting the
   calculation of the dynamical Green’s function. ``"Normal"`` for
   calculating that function from scratch, ``"NoIteration"`` for
   calculating that function with the same iteration in the previous run
   (In this case, the Hamiltonian-vector product is not performed.
   Although the numerical cost is very small, the convergence is not
   guaranteed), ``"Restart_out"`` for calculating that function from
   scratch and writing the restart-file at the end, ``"Restart_in"`` for
   starting the calculation with the previously written restart-file,
   ``"Restart"`` for ``"Restart_out"`` + ``"Restart_in"``.

   The scheme for the spectrum calculation is specified by using the
   parameter ``method``. If ``method="CG"`` is chosen, the shifted
   bi-conjugate gradient method [#]_ together
   with the seed-switch technique
   [#]_ is employed with the
   help of the :math:`K\omega` library [#]_.

* ``SpectrumType``

   **Type :** String (choose from ``"SzSz"``, ``"S+S-"``, ``"Density"``,
   ``"up"``, ``"down"``. ``"SzSz"`` as default.)

   **Description :** The type of the dynamical Green’s function to be
   computed is specified. ``"SzSz"`` for
   :math:`\langle {\hat S}_{z q} {\hat S}_{z q}\rangle`, ``"S+S-"`` for
   :math:`\langle {\hat S}^{+}_{q} {\hat S}^{-}_{q}\rangle`,
   ``"Density"`` for :math:`\langle {\hat n}_{q} {\hat n}_{q}\rangle`,
   ``"up"`` for
   :math:`\langle {\hat c}^{\dagger}_{q \uparrow} {\hat c}_{q \uparrow}\rangle`,
   ``"down"`` for
   :math:`\langle {\hat c}^{\dagger}_{q \downarrow} {\hat c}_{q \downarrow}\rangle`.

*  ``SpectrumQW``, ``SpectrumQL``

   **Type :** Double (default value: ``0.0``)

   **Description :** The wave number (Fractional coordinate) of the
   dynamical Green’s function is specified. The reciprocal lattice
   vector is computed from the direct lattice vector shown in
   :numref:`fig_chap04_1_lattice` , :numref:`fig_chap04_1_honeycomb` ,
   :numref:`fig_ladder` , :numref:`fig_kagome` .

*  ``OmegaMin``

   **Type :** Double (``-LargeValue`` times the number of sites as
   default.)

   **Description :** The lower limit of the real part of the frequency.

*  ``OmegaMax``

   **Type :** Double (``LargeValue`` times the number of sites as
   default.)

   **Description :** The upper limit of the real part of the frequency.

*  ``OmegaIm``

   **Type :** Double (``0.01*LargeValue`` as a default.)

   **Description :** The imaginary part of the frequency.

*  ``NOmega``

   **Type :** Positive integer (``200`` as a default.)

   **Description :** The number of frequencies.
   
.. [#] \A. Frommer, Computing **70**, 87{109 (2003).
.. [#] \S. Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
.. [#] https://github.com/issp-center-dev/Komega.

.. raw:: latex

   \newpage
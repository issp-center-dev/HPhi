.. highlight:: none

Parameters for the dynamical Green’s function
---------------------------------------------

* ``CalcSpec``

   **Type :** String(choose from ``"None"``, ``"Normal"``,
   ``"NoIteration"``, ``"Restart_out"``, ``"Restart_in"``,
   ``"Restart"``, ``"None"`` as default.)

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
   ``"up"``, ``"down"``, ``"SzSz_R"``, ``"S+S-_R"``, ``"Density_R"``, ``"up_R"``,
   ``"down_R"``. ``"SzSz"`` as default.)

   **Description :** The type of the dynamical Green’s function to be
   computed is specified.
   The following values are used For the correlation function in the reciplocal space:
   ``"SzSz"`` for
   :math:`\langle {S}^z_{-\bf q} {S}^z_{\bf q}\rangle`, ``"S+S-"`` for
   :math:`\langle {S}^{+}_{-\bf q} {S}^{-}_{\bf q}\rangle`,
   ``"Density"`` for :math:`\langle {n}_{-\bf q} {n}_{\bf q}\rangle`,
   ``"up"`` for
   :math:`\langle {c}^{\dagger}_{{\bf q} \uparrow} {c}_{{\bf q} \uparrow}\rangle`,
   ``"down"`` for
   :math:`\langle {c}^{\dagger}_{{\bf q} \downarrow} {c}_{{\bf q} \downarrow}\rangle`.
   For the real space, the following values are used:
   ``"SzSz_R"`` for :math:`\langle {\hat S}_{z R} {\hat S}_{z 0}\rangle`,
   ``"S+S-_R"`` for :math:`\langle {\hat S}^{+}_{R} {\hat S}^{-}_{0}\rangle`,
   ``"Density_R"`` for :math:`\langle {\hat n}_{R} {\hat n}_{0}\rangle`,
   ``"up_R"`` for :math:`\langle {\hat c}^{\dagger}_{R \uparrow} {\hat c}_{0 \uparrow}\rangle`, and
   ``"down_R"`` for :math:`\langle {\hat c}^{\dagger}_{R \downarrow} {\hat c}_{0 \downarrow}\rangle`.
   Here :math:`R` spans all site index.
   See :ref:`Fourier-Transformation utility <fourier>` to compute dynamical correlation function in the reciplocal space with the Fourier transformation.

*  ``SpectrumQW``, ``SpectrumQL``

   **Type :** Double (default value: ``0.0``)

   **Description :** The wave number (Fractional coordinate) of the
   dynamical Green’s function is specified. The reciprocal lattice
   vector is computed from the direct lattice vector shown in
   :numref:`fig_chap04_1_lattice` , :numref:`fig_chap04_1_honeycomb` ,
   :numref:`fig_kagome` , :numref:`fig_ladder` .

*  ``OmegaOrg``

   **Type :** Double (``0.0`` as default.)

   **Description :** The origin value of the frequency.

*  ``OmegaMin``

   **Type :** Double (``-LargeValue`` times the number of sites as
   default.)

   **Description :** The lower limit of the real part of the frequency measured from ``OmegaOrg``.

*  ``OmegaMax``

   **Type :** Double (``LargeValue`` times the number of sites as
   default.)

   **Description :** The upper limit of the real part of the frequency measured from ``OmegaOrg``.


*  ``OmegaIm``

   **Type :** Double (``0.01*LargeValue`` as a default.)

   **Description :** The imaginary part of the frequency.

*  ``NOmega``

   **Type :** Positive integer (``200`` as a default.)

   **Description :** The number of frequencies.
   
.. [#] \A. Frommer, Computing **70**, 87-109 (2003).
.. [#] \S. Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
.. [#] https://github.com/issp-center-dev/Komega.

.. raw:: latex

   \newpage

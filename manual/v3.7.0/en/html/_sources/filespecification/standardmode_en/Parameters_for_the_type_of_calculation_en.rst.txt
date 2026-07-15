.. highlight:: none

Parameters for the type of calculation
--------------------------------------

*  ``model``

   **Type :** String (choose from ``"Fermion Hubbard"``, ``"Spin"``,
   ``"Kondo Lattice"``, ``"Fermion HubbardGC"``, ``"SpinGC"``,
   ``"Kondo LatticeGC"``, ``"SpinGCCMA"``) [#]_

   **Description :** The target model is specified with this parameter;
   the expressions above denote the canonical ensemble of the Fermion in
   the Hubbard model

   .. math::
      :label: fml4_1_hubbard

      \mathcal H = -\mu \sum_{i \sigma} c^\dagger_{i \sigma} c_{i \sigma} 
      - \sum_{i \neq j, \sigma} t_{i j} c^\dagger_{i \sigma} c_{j \sigma} 
      + \sum_{i} U n_{i \uparrow} n_{i \downarrow}
      + \sum_{i \neq j} V_{i j} n_{i} n_{j},

   the canonical ensemble in the Spin
   model(\ :math:`\{\alpha, \beta\}=\{x, y, z\}`)

   .. math::
      :label: fml4_1_spin

      \mathcal H &= -h \sum_{i} S^z_{i} - \Gamma \sum_{i} S^x_{i} + D \sum_{i} S^z_{i} S^z_{i}
      \nonumber \\
      &+ \sum_{i j, \alpha}J_{i j \alpha} S^\alpha_{i} S^\alpha_{j}+ \sum_{i j, \alpha \neq \beta} J_{i j \alpha \beta} S_{i}^\alpha S_{j}^\beta,

   the canonical ensemble in the Kondo lattice model

   .. math::
      :label: fml4_1_kondo

      \mathcal H = - \mu \sum_{i \sigma} c^\dagger_{i \sigma} c_{i \sigma} 
      - t \sum_{\langle i j \rangle \sigma} c^\dagger_{i \sigma} c_{j \sigma} 
      + \frac{J}{2} \sum_{i} \left\{
      S_{i}^{+} c_{i \downarrow}^\dagger c_{i \uparrow}
      + S_{i}^{-} c_{i \uparrow}^\dagger c_{i \downarrow}
      + S_{i}^z (n_{i \uparrow} - n_{i \downarrow})
      \right\},

   the grand canonical ensemble of the Fermion in the Hubbard model
   [Eqn. :eq:`fml4_1_hubbard` ], the grand canonical
   ensemble in the Spin model [Eqn. :eq:`fml4_1_spin` ],
   and the grand canonical ensemble in the Kondo lattice model 
   [Eqn. :eq:`fml4_1_kondo` ], respectively.

   When ``model="SpinGCCMA"``, by using a more efficient algorithm [#]_,
   :math:`{\mathcal H}\Phi` calculates a system that is the same as ``"SpinGC"``.
   However, supported models and MPI processes are highly limited. See
   ``"Lattice"`` section.

-  ``method``

   **Type :** String (choose from ``"Lanczos"``, ``"TPQ"``,
   ``"Full Diag"``, ``"CG"``, ``"Time Evolution"``)

   **Description :** The calculation type is specified with this
   parameter; the above expressions above denote the single eigenstate
   calculation by using the Lanczos method, at the finite-temperature by
   using the thermally pure quantum state, the full diagonalization
   method, the multiple eigenstates calculation by using the LOBCG
   method [#]_ [#]_ ,
   and the simulation of real-time evolution, respectively.

   The scheme employed for the spectrum calculation is also specified
   with this parameter. If ``"CG"`` is chosen, the shifted bi-conjugate
   gradient method [#]_ together with the
   seed-switch technique [#]_
   is employed with the help of the :math:`K\omega` library
   [#]_ .

*  ``lattice``

   **Type :** String (choose from ``"Chain Lattice"``,
   ``"Square Lattice"``, ``"Triangular Lattice"``,
   ``"Honeycomb Lattice"``, ``"Kagome"``, ``"Ladder"``)

   **Description :** The lattice shape is specified with this parameter;
   the expressions above denote the one-dimensional chain lattice ( :numref:`fig_chap04_1_lattice` (a)),
   the two-dimensional square lattice ( :numref:`fig_chap04_1_lattice` (b)),
   the two-dimensional triangular lattice ( :numref:`fig_chap04_1_lattice` (c)),
   the two-dimensional anisotropic honeycomb lattice ( :numref:`fig_chap04_1_honeycomb` ),
   the Kagome Lattice( :numref:`fig_kagome` ),
   and the ladder lattice ( :numref:`fig_ladder` ) respectively.

   In ``method="SpinGCCMA"``, only ``"Chain Lattice"``,
   ``"Honeycomb Lattice"``, ``"Kagome"``, and ``"Ladder"`` are
   supported. The limits of :math:`L`, :math:`W`, and the number of MPI
   processes (:math:`N_{\rm proc}`) are as follows:

   *  ``"Chain Lattice"``

      :math:`L = 8n` (where :math:`n` is an integer number under the
      condition :math:`n\geq1`), :math:`N_{\rm proc} \leq 2(L=8)`,
      :math:`N_{\rm proc} \leq 2^{L/2-2}(L>8)`.

   *  ``"Honeycomb Lattice"``

      :math:`W=3, L \geq 2`, :math:`N_{\rm proc} \leq 2(L=2)`,
      :math:`N_{\rm proc} \leq 64(L>2)`.

   *  ``"Kagome"``

      :math:`W=3, L \geq 2`, :math:`N_{\rm proc} \leq 1(L=2)`,
      :math:`N_{\rm proc} \leq 512(L>2)`.

   *  ``"Ladder"``

      :math:`W=2, L = 2n` (where :math:`n` is an integer number under
      the condition :math:`n\geq4`), :math:`N_{\rm proc} \leq 2^{L-4}`.

.. [#] \GC=Grand Canonical
.. [#] \Y. Yamaji *et. al.*, manuscript in preparation.
.. [#] A.V.Knyazev, SIAM Journal on Scientific Computing **23**, 517 (2001).
.. [#] S.Yamada, T.Imamura, M.Machida, The Japan Society for Computational Engineering and Science **2006**, 20060027 (2006).
.. [#] A.Frommer, Computing **70**, 87-109 (2003).
.. [#] S.Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
.. [#] https://github.com/issp-center-dev/Komega.


.. raw:: latex

   \newpage

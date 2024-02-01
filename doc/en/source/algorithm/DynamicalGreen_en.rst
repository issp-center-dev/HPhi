.. highlight:: none

Dynamical Green's function
--------------------------

Using :math:`{\mathcal H}\Phi`, we can calculate a dynamical Green's function

.. math:: G_n^{O_l,O_r}(z) = \langle \Phi_n | \hat{O}_l (z + E_n - \hat{\cal H})^{-1} \hat{O}_r| \Phi_n \rangle

where :math:`\hat{O}_{l,r}` is a single exciation operator

.. math:: \sum_{i, \sigma_1} A_{i \sigma_1} c_{i \sigma_1} \quad \textrm{or} \quad \sum_{i, \sigma_1} A_{i \sigma_1} c_{i\sigma_1}^{\dagger}

or a pair-exciation operator

.. math:: \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i \sigma_1}c_{j \sigma_2}^{\dagger} \quad \textrm{or} \quad
          \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i\sigma_1}^{\dagger}c_{j\sigma_2}.

For example, to compute the dynamical spin susceptibility, we use pair excitation operators

.. math:: \hat{O}_r = \hat{S}_{\textbf{R}=\textbf{0}}^z = \frac{1}{2} (c_{\textbf{0}\uparrow}^{\dagger}c_{\textbf{0}\uparrow}-c_{\textbf{0}\downarrow}^{\dagger}c_{\textbf{0}\downarrow})
    \\
    \hat{O}_l = \hat{S}_{\textbf{R}}^z = \frac{1}{2} (c_{\textbf{R}\uparrow}^{\dagger}c_{\textbf{R}\uparrow}-c_{\textbf{R}\downarrow}^{\dagger}c_{\textbf{R}\downarrow})

to generate :math:`G_n^{O_l,O_r}(z)\equiv G_n^{\textbf{R}}(z)`,
then perform the Fourier transformation

.. math:: G_n^{\textbf{k}}(z) \equiv \sum_{\textbf{R}} \exp(i\textbf{k}\cdot\textbf{R}) G_n^{\textbf{R}}(z)

as a postprocess.

Three modes are implemented in :math:`{\cal H}\Phi`:
The continued fraction expansion method by using Lanczos method [1]_,
the shifted Krylov method [2]_, and
the Lehmann representation with the full diagonallization

.. math:: G_n^{O_l,O_r}(z) = \sum_{m} \frac{\langle \Phi_n | \hat{O}_l | \Phi_m \rangle \langle \Phi_m |\hat{O}_r| \Phi_n \rangle}{z + E_n - E_m}.

See the reference for the details of each algorithm.

.. [1] \E. Dagotto, Rev. Mod. Phys. **66**, 763-840 (1994).
.. [2] \S.Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).

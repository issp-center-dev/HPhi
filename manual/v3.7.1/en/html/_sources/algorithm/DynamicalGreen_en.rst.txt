.. highlight:: none

Dynamical Green’s function
--------------------------

Using :math:`{\mathcal H}\Phi`, we can calculate a dynamical Green’s
function

.. math:: I(z) = \langle \Phi ' | \frac{1}{ {\mathcal H}- z\hat{I} } | \Phi '\rangle,

where :math:`|\Phi ' \rangle  = \hat{O} | \Phi _0 \rangle` is an
excited state and :math:`\hat{O}` is an excitation operator defined as a
single excitation operator

.. math:: \sum_{i, \sigma_1} A_{i \sigma_1} c_{i \sigma_1} (c_{i\sigma_1}^{\dagger})

or a pair excitation operator

.. math:: \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i \sigma_1}c_{j \sigma_2}^{\dagger} (c_{i\sigma_1}^{\dagger}c_{j\sigma_2}).

For example, the dynamical spin susceptibilities can be calculated by
defining :math:`\hat{O}` as

.. math:: \hat{O} = \hat{S}({\bf k}) = \sum_{j}\hat{S}_j^z e^{i  {\bf k} \cdot \bf {r}_j} = \sum_{j}\frac{1}{2} (c_{j\uparrow}^{\dagger}c_{j\uparrow}-c_{j\downarrow}^{\dagger}c_{j\downarrow})e^{i  {\bf k} \cdot \bf {r}_j}.

There are two modes implemented in :math:`{\cal H}\Phi`. One is the
continued fraction expansion method by using Lanczos method
 [#]_ and the other is the shifted Krylov
method [#]_ . See the reference
for the details of each algorithm.

.. [#] \E. Dagotto, Rev. Mod. Phys. **66**, 763-840 (1994).
.. [#] \S.Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
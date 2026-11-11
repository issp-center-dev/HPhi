Overview
========

In this document, we introduce how we compute downfolded models
with mVMC or :math:`{\mathcal H}\Phi` in conjunction to
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_.
In RESPACK, the screened direct integrals :math:`U_{mn}({\bf R},\omega)` and the screened exchanged integrals :math:`J_{mn}({\bf R},\omega)` are given as follows:

.. math::
   \begin{aligned}
   U_{mn}({\bf R},\omega)&=&\int_V d{\bf r} \int_V  d{\bf r'}
   w_{m{\bf 0}}^*({\bf r}) w_{m{\bf 0}}({\bf r}) 
   W({\bf r,r'},\omega)
   w_{n{\bf R}}^*({\bf r'}) w_{n{\bf R}}({\bf r'}),\nonumber\\
   J_{mn}({\bf R},\omega)&=&\int_V  d{\bf r}\int_V d{\bf r'}
   w_{m{\bf 0}}^*({\bf r}) w_{n{\bf R}}({\bf r}) 
   W({\bf r,r'},\omega) 
   w_{n{\bf R}}^*({\bf r'}) w_{m{\bf 0}}({\bf r'}). 
   \end{aligned}

Here, :math:`V` is the volume of the crystal, :math:`w_ {i {\bf R}}({\bf r})` is the :math:`i` -th wannier function at :math:`\bf R` -th cell, :math:`W({\bf r,r'}, \omega)` is the screened Coulomb interactions, respectively. In the following, the components at :math:`\omega=0` are only considered. Then, the Hamiltonian of the two-body interactions are given as follows:

.. math::
   \begin{aligned}
   {\cal H}_{\rm int} &= \frac{1}{2}\sum_{\sigma\rho }\sum_{ij}\sum_{nm} \Bigl[ U_{mn}({\bf R}_{ij})c_{im, \sigma}^{\dagger}c_{jn, \rho}^{\dagger}c_{jn, \rho}c_{im, \sigma}\nonumber\\
   &+ J_{mn}({\bf R}_{ij})(c_{im, \sigma}^{\dagger}c_{jn,\rho}^{\dagger}c_{im,\rho}c_{jn,\sigma} + c_{im, \sigma}^{\dagger}c_{im,\rho}^{\dagger}c_{jn,\rho}c_{jn,\sigma}  )\Bigr],
   \end{aligned}

where :math:`{\bf R}_{ij} \equiv {\bf R}_i-{\bf R}_j` . Since mVMC and :math:`{\mathcal H}\Phi` cannot directly treat the following type of interactions :math:`{c_{i, \sigma}^{\dagger}c_{j, \rho}^{\dagger}c_{k, \rho'}c_{l, \sigma'}}` , the Hamiltonian must be rewritten as follows:

.. math::
   \begin{aligned}
   {\cal H}_{\rm int} &= \sum_{i,m} U_{mm}({\bf 0})n_{im,\uparrow} n_{im, \downarrow} +\sum_{(i,m)<(j,n)}U_{mn}({\bf R}_{ij})n_{im}n_{jn}\nonumber\\
   & - \sum_{(i,m)<(j,n)}J_{mn}({\bf R}_{ij})(n_{im, \uparrow}n_{jn,\uparrow}+n_{im, \downarrow}n_{jn,\downarrow}) \nonumber\\
   & + \sum_{(i,m)<(j,n)}J_{mn}({\bf R}_{ij})(c_{im, \uparrow}^{\dagger}c_{jn,\downarrow}^{\dagger}c_{im,\downarrow}c_{jn,\uparrow}+{\rm h.c.}) \nonumber\\
   & + \sum_{(i,m)<(j,n)}J_{mn}({\bf R}_{ij}) (c_{im, \uparrow}^{\dagger}c_{im,\downarrow}^{\dagger}c_{jn,\downarrow}c_{jn,\uparrow} + {\rm h.c.} ).
   \end{aligned}

The lattice model is defined by the following Hamiltonian:
   
.. math::

   \begin{aligned}
   {\cal H} &=
   \sum_{m,n, i, j,\sigma}
   \left[t_{mn}({\bf R}_{ij}) - t_{mn}^{\rm DC}({\bf R}_{ij})\right] c_{im \sigma}^{\dagger} c_{jn \sigma}
   + {\cal H}_{int},
   \end{aligned}

where  :math:`t_{mn}^{\rm DC}({\bf R}_{ij})` is the one-body correction term given as:

.. math::

   \begin{aligned}
   t_{mm}^{\rm DC}({\bf 0}) &\equiv \alpha U_{mm}({\bf 0}) D_{mm}({\bf 0})
   + \sum_{({\bf R}, n) \neq ({\bf 0}, m)} U_{m n} ({\bf R})D_{nn}({\bf 0})\\
   & - (1-\alpha) \sum_{({\bf R}, n) \neq ({\bf 0}, 0)} J_{m n}({\bf R}) D_{nn}({\bf R}),\\
   t_{mn}^{\rm DC}({\bf R}_{ij}) &\equiv \frac{1}{2} J_{mn}({\bf R}_{ij}) \left(D_{nm}({\bf R}_{ji}) + 2 {\rm Re} [D_{nm}({\bf R}_{ji})]\right)\\
   &-\frac{1}{2}  U_{mn}({\bf R}_{ij}) D_{nm}({\bf R}_{ji}),
   \quad ({\bf R}_{ij}, m) \neq ({\bf 0}, n),
   \\
   D_{mn}({\bf R}_{ij}) &\equiv \sum_{\sigma}
   \left\langle c_{im \sigma}^{\dagger} c_{jn \sigma}\right\rangle_{\rm KS},
   \end{aligned}

Here, :math:`t_{mm}^{\rm DC}({\bf 0})` is the term to correct the chemical potntial, :math:`t_{mn}^{\rm DC}({\bf R}_{ij})` is term to correct transfer integrals
. These terms are introduced to avoid double counting in analyzing the lattice model. To adopt theses corrections or not can be selected by the option ``doublecounting`` in the input file. The strength of :math:`U_{Rij}` and :math:`J_{Rij}` can be controled by multiplying tuning parameters :math:`\lambda_U, \lambda_J`. For details, see ``Input parameters for Standard mode``.

Prerequisite
------------

We compute the Kohn-Sham orbitals with
`QuantumESPRESSO <http://www.quantum-espresso.org/>`_
or
`xTAPP <http://xtapp.cp.is.s.u-tokyo.ac.jp/>`_,
and obtain the Wannier function, the dielectric function,
the effective interaction with RESPACK,
and simulate quantum lattice models with
mVMC or :math:`{\mathcal H}\Phi`.
Therefore, these programs must be available in our machine.

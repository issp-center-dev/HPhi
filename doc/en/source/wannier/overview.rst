Overview
========

In this document, we introduce how we compute downfolded models
with mVMC or :math:`{\mathcal H}\Phi` in conjunction to
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_.

.. math::

   \begin{aligned}
   {\cal H} &=
   \sum_{i, j, \alpha, \beta, \sigma}
   t_{i \alpha j \beta} c_{i \alpha \sigma}^{\dagger} c_{j \beta \sigma}
   \nonumber \\
   &+ \sum_{i, \alpha}
   U_{i \alpha i \alpha} n_{i \alpha \uparrow} n_{j \alpha \downarrow}
   + \sum_{(i, \alpha) \lt (j, \beta)}
   U_{i \alpha j \beta} n_{i \alpha} n_{j \beta}
   - \sum_{(i, \alpha) \lt (j, \beta)}
   J_{i \alpha j \beta} (n_{i \alpha \uparrow} n_{j \beta \uparrow}
   + n_{i \alpha \downarrow} n_{j \beta \downarrow})
   \nonumber \\
   &+ \sum_{(i, \alpha) \lt (j, \beta)}
   J_{i \alpha j \beta} (
   c_{i \alpha \uparrow}^{\dagger} c_{j \beta \downarrow}^{\dagger}
   c_{i \alpha \downarrow} c_{j \beta \uparrow} +
   c_{j \beta \uparrow}^{\dagger} c_{i \alpha \downarrow}^{\dagger}
   c_{j \beta \downarrow} c_{j \alpha \uparrow} )
   \nonumber \\
   &+ \sum_{(i, \alpha) \lt (j, \beta)}
   J_{i \alpha j \beta} (
   c_{i \alpha \uparrow}^{\dagger} c_{i \alpha \downarrow}^{\dagger}
   c_{j \beta \downarrow} c_{j \beta \uparrow} +
   c_{j \beta \uparrow}^{\dagger} c_{j \beta \downarrow}^{\dagger}
   c_{i \alpha \downarrow} c_{i \alpha \uparrow} )
   \end{aligned}

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

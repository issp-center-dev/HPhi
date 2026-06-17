.. highlight:: none

.. _Sec:sec_bogoliubov_rep:

Bogoliubov representation
=========================


In the spin system,
the spin indices in the input files of ``transfer``, ``InterAll``,
and correlation functions are specified as those of the Bogoliubov representation.
The spin operators are written by using creation\/annihilation operators:

.. math::

  S_{i z} &= \sum_{\sigma = -S}^{S} \sigma c_{i \sigma}^\dagger c_{i \sigma}
  \\
  S_{i}^+ &= \sum_{\sigma = -S}^{S-1} 
  \sqrt{S(S+1) - \sigma(\sigma+1)} 
  c_{i \sigma+1}^\dagger c_{i \sigma}
  \\
  S_{i}^- &= \sum_{\sigma = -S}^{S-1} 
  \sqrt{S(S+1) - \sigma(\sigma+1)} 
  c_{i \sigma}^\dagger c_{i \sigma+1}

In HPhi, the index of the highest-:math:`\sigma` state is 0.

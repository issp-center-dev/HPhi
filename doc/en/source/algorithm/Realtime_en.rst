.. highlight:: none

Real time evolution method
==========================

In :math:`{\mathcal H}\Phi`, real time evolution calculation is done by using the following relation

.. math::

 |\Phi (t_n)\rangle = \exp (-i {\cal H}  \Delta t_n)|\Phi (t_{n-1})\rangle,

where :math:`|\Phi(t_0)\rangle` is an initial wave function and :math:`t_n = \sum_{j=1}^n  \Delta t_j`.
In calculation, we approximate :math:`\exp (-i {\cal H}  \Delta t_n)` as

.. math::

 \exp (-i {\cal H}  \Delta t_n) =\sum_{l=0}^m \frac{1}{l!}(-i {\cal H}  \Delta t_n)^l .

Here, the cut-off integer :math:`m` can be set by `ExpandCoef` in `ModPara`.
We can judge whether the expansion order is enough or not by checking the norm conservation :math:`\langle \Phi (t_n)|\Phi (t_n)\rangle=1` and energy conservation :math:`\langle \Phi (t_n)|\hat{\cal H}|\Phi (t_n)\rangle=E`.

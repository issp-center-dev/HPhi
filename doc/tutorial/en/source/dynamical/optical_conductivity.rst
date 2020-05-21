Hubbard chain (optical conductivity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here, we calculate the optical conductivity
for the one-dimensional Hubbard model.

The optical conductivity :math:`\sigma(\omega)` can be calculated from
the current-current correlation 
:math:`I(\omega,\eta)`, which is defined as

.. math::

 j_{x}={i}\sum_{i,\sigma}(c_{{\bf r}_{i}+{\bf e}_{x},\sigma}^{\dagger}c_{{\bf r}_{i},\sigma}-c_{{\bf r}_{i},\sigma}^{\dagger}c_{{\bf r}_{i}+{\bf e}_{x},\sigma}), 

 I(\omega,\eta)={\rm Im}\Big[\langle 0|j_{x}[H-(\omega-E_{0}-{i}\eta)I]^{-1}j_{x}|0\rangle\Big],

where :math:`{\bf e}_{x}` is the unit translational vector
in the x direction.
From this
the regular part of the optical conductivity 
is defined as

.. math::

 \sigma_{\rm reg}(\omega)=\frac{I(\omega,\eta)+I(-\omega,-\eta)}{\omega N_{s}},

where :math:`N_{\rm s}` is the number of sites.

An input file (stan.in) for 6-site Hubbard model is as follows::

 model = "Hubbard" 
 method = "CG" 
 lattice = "chain" 
 L = 6
 t = 1
 U = 10
 2Sz = 0
 nelec = 6
 exct = 1
 EigenVecIO  = "out"


Scripts for calculating the optical conductivity are 
available at **HPhi/tool/ForOpticalConductivity**.
  
By performing the all-in-one script (All.sh),  ::

 sh ./All.sh

you can obtain **optical.dat**.

A way for plotting **optical.dat** is as follows  ::

 plot "optical.dat" u 1:(-($4+$8)/$1) w l


Tutorial for finite-temperature calculations
==============================
Heisenberg chain (finite temperatures)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here, we study the finite-temperature
properties of spin 1/2 Heisenberg model on the chain.

.. math::

 H = J \sum_{\langle i,j\rangle}{\bf S}_{i}\cdot{\bf S}_{j}

The input file (stan.in) for 12-site Heisenberg model is as follows::

 L       = 12
 model   = "Spin" 
 method  = "FullDiag" 
 lattice = "chain"
 J = 1
 2Sz = 0
 2S  = 1

You can execute HPhi as follows ::

 HPhi -s stan.in

Full diagonalization
"""""""""""""""""""""""""""""""
After executing the full diagonalization,
all the eigen energies are output in **output/Eiegenvalue.dat**.
By using the python script **(Git/HPhi/tool/FiniteT/Finite.py)**, 
you can obtain the temperature dependence of the energy and the specific heat.

You can execute **Finite.py** as follows ::

 python3 Finite.py

Then, you can obtain **FullDiag.dat** as follows ::

     0.000100  -5.3873909174000003    0.0000000000000000   
     0.000150  -5.3873909174000003    0.0000000000000000   
     0.000200  -5.3873909174000003    0.0000000000000000   
     0.000250  -5.3873909174000003    0.0000000000000000   
     0.000300  -5.3873909174000003    0.0000000000000000   

The 1st row represents temperature, 2nd row represents the energy, and
the 3rd row represents the specific heat defined 
by :math:`C=(\langle E^2 \rangle-\langle E \rangle^2)/T^2`.

TPQ method
"""""""""""""""""""""""""""""""
By selecting method as "TPQ",
you can perform the finite-temperature calculations using the TPQ method.

The input file (stan.in) for 12-site Heisenberg model is as follows::

 L       = 12
 model   = "Spin" 
 method  = "TPQ" 
 lattice = "chain"
 J = 1
 2Sz = 0
 2S  = 1

After performing the TPQ calculations,
all the eigen energies are output in **output/SS_rand*.dat**.
By using the python script **(Git/HPhi/tool/FiniteT/AveSSrand.py)**, 
you can obtain the temperature dependence of 
physical quantities such as the energy and the specific heat.

Then, you can obtain **ave_TPQ.dat** as follows ::

 # temperature T    err of T  Energy E    error of E  specific heat C  err of C   
 30.17747           0.02022   -0.35489    0.04044     0.00264          8.1126-05
 15.10875           0.01016   -0.43500    0.04065     0.01068          0.0003182
 10.08598           0.00681   -0.51592    0.04086     0.02423          0.0007030
 7.574693           0.00513   -0.59754    0.04106     0.04335          0.0012311
 6.067978           0.00412   -0.67978    0.04123     0.06808          0.0019053

Using gnuplot, you can directly compare the two results :: 

  plot "FullDiag.dat" u 1:2 w l,"ave_TPQ.dat" u 1:3:4 w e

You can see the following output image.

.. image:: ../../figs/finiteT.*
   :align: center

Kitaev cluster (finite temperatures)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

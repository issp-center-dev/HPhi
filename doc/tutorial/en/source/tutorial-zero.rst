Tutorial for calculations at zero temperature
==============================
Spin 1/2 Dimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following spin 1/2 dimer model (2-site Heisenberg model).

.. math::

 H = J {\bf S}_{0}\cdot{\bf S}_{1}

The input file (stan.in) is as follows::

 L=2
 model = "Spin" 
 method = "FullDiag" 
 lattice = "chain"
 J = 0.5
 2Sz = 0
 2S  = 1

You can execute HPhi as follows ::

 HPhi -s stan.in

Check the energy
"""""""""""""""""""""""""""""""
Please check whether the energies are given as follows.

 :math:`E_{\rm min}=-3/4` (singlet) 

 :math:`E_{\rm max}=1/4` (triplet) 

Check S dependence
"""""""""""""""""""""""""""""""
By changing 2S=1 in stan.in, you can treat spin-S dimer (eg. 2S=2 means S=1).
Please check whether the energies are given as follows.

 :math:`E_{\rm min}=-S(S+1)` 

 :math:`E_{\rm max}=S^2` 

Add magnetic field H
"""""""""""""""""""""""""""""""
By adding  H in stan.in, selecting model as "SpinGC", and deleting "2Sz=0" 
you can examine the 
effects of the magnetic field in the dimer model.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "FullDiag" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2

Please check whether the ground state becomes polarized state (Sz=1).


Try to use Lanczos method
"""""""""""""""""""""""""""""""
**This is just a pedagogical example.**

By selecting method as "Lanczos"
you can perform the Lanczos calculations.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "Lanczos" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2

This calculation will **fail**! 
Please think why the Lanczos method fails for the dimer. 


Try to use LOBCG method
"""""""""""""""""""""""""""""""
LOBCG is locally optimal block conjugate gradient method.
By selecting method as "CG",
you can perform the LOBCG calculations.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "CG" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2


In contrast to the Lanczos method, 
this calculation will work well ! 
Please think why the CG method works well for the dimer. 

Please also check the excited states can
be correctly obtained by using LOBCG method.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "CG" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2
 exct = 4

Here, exct represents the number of excited states, which are
obtained by the LOBCG method.

Hubbard Dimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following the Hubbard dimer model.

.. math::

 H = -t \sum_{\sigma}(c_{0\sigma}^{\dagger}c_{1\sigma}+{\rm H.c.})
   +U(n_{0\uparrow}n_{0\downarrow}+n_{1\uparrow}n_{1\downarrow})

The input file (stan.in) is as follows::

 model = "Hubbard" 
 method = "FullDiag" 
 lattice = "chain" 
 L=2
 t = -0.5 
 U = 4
 2Sz = 0
 nelec = 2

You can execute HPhi as follows ::

 HPhi -s stan.in

Check the energy
"""""""""""""""""""""""""""""""
For the Hubbard dimer at half filling with total Sz=0, 
energies are given as follows:

 :math:`E=0,U,\frac{U}{2}\times(1\pm\sqrt{(1+(4t/U)^2)})` 

For example, by taking :math:`U=4,t=-1`, the 
energies are  given as follows:

 :math:`E=-0.828427, 0, 4, 4.828427` 

It is note that simple mathematical calculations 
can be done using:: 

 bc -l 

on the terminal.

Try to use LOBCG method
"""""""""""""""""""""""""""""""
The input file (stan.in) is as follows::

 model = "Hubbard" 
 method = "CG" 
 lattice = "chain" 
 L=2
 t = -0.5 
 U = 4
 2Sz = 0
 nelec = 2
 exct = 4

Please check whether LOBCG method correctly 
reproduces the energies including the excited states.

Hubbard Trimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following the Hubbard trimer model 
(Hubbard model on a triangle).

.. math::

 H = -t \sum_{\sigma}(c_{0\sigma}^{\dagger}c_{1\sigma}+c_{1\sigma}^{\dagger}c_{2\sigma}
   +c_{2\sigma}^{\dagger}c_{0\sigma}+{\rm H.c.})
   +U\sum_{i}(n_{i\uparrow}n_{i\downarrow})

The input file (stan.in) is as follows::

 model = "Hubbard" 
 method = "FullDiag" 
 lattice = "chain" 
 L=2
 t = -1
 U = 4
 2Sz = 0
 nelec = 2

Note that the filling is not half filling.

You can execute HPhi as follows ::

 HPhi -s stan.in

Ferromagnetic ground state
"""""""""""""""""""""""""""""""
For the Hubbard model on a triangle with one hole, 
it is known that the **perfect ferromagnetism** becomes ground state.
Please check that. 

If you want know the mechanism of the 
ferromagnetism, please see 
**Hal Tasaki, Kotai Butsuri, Vol. 31, 173 (1996)**.
This is one of the simplest example of the 
Nagaoka's ferromagnetism.


Effects of transfer integrals
"""""""""""""""""""""""""""""""
Please the effects of the sign of
the transfer integrals. 
**For example, what happens if you take t = 1 ?**

Another interesting example is by changing 
the transfer integrals between site 0 and site 2.
Following an example of the **trans.def** ::

  ======================== 
  NTransfer      12  
  ======================== 
  ========i_j_s_tijs====== 
  ======================== 
    1     0     0     0         -1.000000000000000         0.000000000000000
    0     0     1     0         -1.000000000000000         0.000000000000000
    1     1     0     1         -1.000000000000000         0.000000000000000
    0     1     1     1         -1.000000000000000         0.000000000000000
    2     0     0     0         -2.000000000000000         0.000000000000000
    0     0     2     0         -2.000000000000000         0.000000000000000
    2     1     0     1         -2.000000000000000         0.000000000000000
    0     1     2     1         -2.000000000000000         0.000000000000000
    2     0     1     0         -1.000000000000000         0.000000000000000
    1     0     2     0         -1.000000000000000         0.000000000000000
    2     1     1     1         -1.000000000000000         0.000000000000000
    1     1     2     1         -1.000000000000000         0.000000000000000

Using this transfer integrals, please examine the
U dependence of the ground state.
Is there phase transition between singlet ground state and
the perfect ferromagnetism ?


Heisenberg chain (zero temperature)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's solve the following spin 1/2 Heisenberg model on the chain.

.. math::

 H = J \sum_{\langle i,j\rangle}{\bf S}_{i}\cdot{\bf S}_{j}

The input file (stan.in) for 16-site Heisenberg model is as follows::

 L       = 16
 model   = "Spin" 
 method  = "CG" 
 lattice = "chain"
 J = 1
 2Sz = 0
 2S  = 1

You can execute HPhi as follows ::

 HPhi -s stan.in

Check the energy
"""""""""""""""""""""""""""""""
Please check whether the energies are given as follows.

.. math::

 E_{0}= -7.142296 

Obtaining the excited state
"""""""""""""""""""""""""""""""
By adding **exct=2**, you can obtain the 2 low-energy states.
Please check the energies.

.. math::
  
 E_{0}= -7.142296

 E_{1}= -6.872107 

Size dependence of the spin gap
"""""""""""""""""""""""""""""""
The spin gap at finite system size is defined
as :math:`\Delta=E_{1}-E_{0}`. For 16-site,
we obtain :math:`\Delta\sim 0.2701`.

Please examine how :math:`\Delta` behaves
as a function of system size L.
(available system size on PC may be L=24)

Haldane gap
"""""""""""""""""""""""""""""""
By performing a similar calculations for S=1 system,
please examine  how :math:`\Delta` behaves
as a function of system size L.
It is known that the finite spin gap exists even
in the thermodynamic limit (:math:`L=\infty`).
This spin gap is often called Haldane gap.

J1-J2 Heisenberg model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Here, we solve the :math:`J_{1}-J_{2}` Heisenberg model on the square lattice, 
which is a canonical example of the frustrated magnets.
Its Hamiltonian is defined as

.. math::

  {\mathcal H}=J_{1}\sum_{\langle i,j\rangle }{\bf S}_{i}\cdot{\bf S}_{j}+J_{2}\sum_{\langle\langle i,j\rangle\rangle }{\bf S}_{i}\cdot{\bf S}_{j},

where :math:`J_{1}` (:math:`J_{2}`) represents the nearest (next-nearest) neighbor interactions.

An input file for treating :math:`J_{1}-J_{2}` Heisenberg model is given as ::

 model = "Spin" 
 method = "CG" 
 lattice = "square" 
 L = 4
 W = 4
 J = 1
 J' = 1
 2S = 1
 2Sz = 0
 exct = 2

Here, J (J') represents :math:`J_{1}` (:math:`J_{2}`).

Calculations of spin structure factors for ground state
"""""""""""""""""""""""""""""""
First, we calculate the spin structure factors, which are defined as

.. math::

  S({\bf q})=\frac{1}{N_{\rm s}}\sum_{i,j} {\bf S}_{i}\cdot{\bf S}_{j}

To calculate :math:`S({\bf q})`, it is necessary to prepare
a proper input file for two-body Green functions.
By using a python script **MakeGreen.py** (HPhi/tool/ForSq/MakeGreen.py),
you can generate **greentwo.def** for calculating :math:`S({\bf q})`.
To use  **MakeGreen.py** , an input file for specifying lattice geometry (**input.txt**) is necessary,
whose form is given as follows ::

 Lx 4
 Ly 4
 Lz 1
 orb_num 1

Here, Lx (orb_num) represents the length of the x direction (number of orbitals).

By using a python script **CalcSq.py** (HPhi/tool/ForSq/CalcSq.py),
you can calculate :math:`S({\bf q})` from **output/zvo_cisajscktalt_eigen0.dat**.

Procedure for calculating and visualizing :math:`S({\bf q})` is given as follows ::
 1. HPhi -sdrt stan.in
 2. pyhton3 MakeGreen.def (input.txt is necessary)
 3. HPhi -e namelist.def 
 4. python3 CalcSq.py (input.txt is necessary)
 5. You can obtain **Sq_eigen0.dat** !!
 6. plot "Sq_eigen0.dat" u 1:2:3 (using gnuplot)

Following the procedure, please see how :math:`S({\bf q})` changes
by changing J'.
As an example, we show :math:`S({\bf q})` for J1/J2=0 and 1 below.

.. image:: ../../figs/J1J2_Sq.*
   :align: center

Calculations of spin structure factors for excited states
"""""""""""""""""""""""""""""""
By changing exct in stan.in, you can obtain several excited states.
For those excited states, by changing **max_num=1** as ,for example, **max_num=4**,
you can obtain  :math:`S({\bf q})` for the excited states.

Please see how :math:`S({\bf q})` in the excited states 
changes by changing J'. For example, please check 
what is the nature of the
first excited state J'=0,0.5, and 1.
 

How to use Expert mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you prepare input files, you can perform calculations for
arbitrary Hamiltonians with any one-body potentials and the two-body interactions.   
By taking spin 1/2 system as an example,
we explain how to prepare input files.
For spin 1/2 system, we prepare simple python scripts (**HPhi/tool/ExpertSpin/MakeDef.py**)that
can generate the input files for general Hamiltonians, which are defined as

.. math::

  {\mathcal H}=\sum_{i,j} J_{i,j}^{\alpha,\beta} {\bf S}_{i}^{\alpha} {\bf S}_{j}^{\beta}.

To use *MakeDef.py*, it is necessary to prepare two input files,
**input.txt** and **pair.txt**. 

In **input.txt**, two parameters **Ns** (number of sites) and **exct** (number of excited states)
are specified.

Below is an example of **input.txt** for 2 site Heisenberg model ::

 Ns 2
 exct 2

In **pair.txt**, you specify the interaction terms in the form

.. math::
  i~~~~~j~~~~~\alpha~~~~~\beta~~~~~J_{i,j}^{\alpha,\beta}


Below is an example of **pair.txt** for 2 site Heisenberg model ::

 0 1 x x 0.5
 0 1 y y 0.5
 0 1 z z 0.5

You can also specify the non-diagonal interaction as ::

 0 1 x x 0.5
 0 1 y y 0.5
 0 1 z z 0.5
 0 1 x y 0.5
 0 1 x z 0.5
 0 1 y z 0.5

Note that interaction terms must be specified for **(x,y), (x,z), (y,z)**
and **(y,x), (z,x), (z,y) cannot be used**.

Exercise
"""""""""""
By changing **pair.txt** and **input.txt**,
you can treat your favorite models.
For example, please try to make input files
for the **Kitaev model** on the honeycomb lattice.
We note that the **Kitaev model** can be used in the Standard mode.

Another example is the **XY model** on the chain.
In the standard model,
you can also treat **XY model** by omitting
"CoulombInter  coulombinter.def"  and
"Hund  hund.def" in namelist.def.


Use eigenvectors
^^^^^^^^^^^^^^^^^^^^^^^^^
In this tutorial, we will study how to read the eigenvectors.
In the standard mode, setting ``EigenVecIO = "Out"`` makes HPhi to write the calculated eigenvectors as ``output/zvo_eigenvec_[index]_rank_[rank].dat``, where ``[index]`` is the index of the states (e.g., the ground state has ``[index] = 0``) and ``[rank]`` is the rank of the process.
In the MPI parallelization with Npara processes, HPhi splits the whole Hilbert space into the Npara blocks and each process treats the one of them.
The file format is described in the `reference manual <http://issp-center-dev.github.io/HPhi/manual/master/en/html/filespecification/outputfiles_en/tmpvec_en.html>`_ .
For example, the following python function reads the vector::

  def read_gs(*, index=0, rank=0):
      import numpy as np
      from os.path import join
      from struct import unpack

      filename = join("output",
                      "zvo_eigenvec_{}_rank_{}.dat".format(index,
                                                           rank))
      with open(filename, "rb") as f:
          f.read(4)
          nelems = unpack("L", f.read(8))[0]
          ret = np.zeros(nelems, dtype=np.complex)
          f.read(16)
          for i in range(nelems):
              re = unpack("d", f.read(8))[0]
              im = unpack("d", f.read(8))[0]
              ret[i] = np.complex(re, im)
          return ret

Exercise
"""""""""""
Check the orthogonality of the eigenvectors calculated by the LOBCG method by calculating the norm and the inner-product of some of the eigenvectors.

Hint: In the standard mode, the ``exct`` keyword controls the number of eigenvectors to be calculated.

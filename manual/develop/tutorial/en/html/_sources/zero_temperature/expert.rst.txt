How to use Expert mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you prepare input files, you can perform calculations for
arbitrary Hamiltonians with any one-body potentials and the two-body interactions.   
By taking spin 1/2 system as an example,
we explain how to prepare input files.
For spin 1/2 system, we prepare simple python scripts (``samples/tutorial_1.6/MakeDef.py``) 
that can generate the input files for general Hamiltonians, which are defined as

.. math::

  {\mathcal H}=\sum_{i,j} J_{i,j}^{\alpha,\beta} {\bf S}_{i}^{\alpha} {\bf S}_{j}^{\beta}.

Note tat ``samples/tutorial_1.6/read.py`` and ``samples/tutorial_1.6/hphi_io.py`` are necessary for **MakeDef.py**.
To use *MakeDef.py*, it is necessary to prepare two input files,
**input.txt** and **pair.txt** (examples of them are available in ``samples/tutorial_1.6`` ). 

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

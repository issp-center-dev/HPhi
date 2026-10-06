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

In **pair.txt**, one specify the interaction terms in the form

.. math::
  i~~~~~j~~~~~\alpha~~~~~\beta~~~~~J_{i,j}^{\alpha,\beta}

For diagonal interactions (:math:`\alpha=\beta`), :math:`J_{i,j}^{\alpha,\alpha} {\bf S}_{i}^{\alpha} {\bf S}_{j}^{\alpha}` is added, 
and :math:`J_{i,j}^{\alpha,\beta} [{\bf S}_{i}^{\alpha} {\bf S}_{j}^{\beta}+{\bf S}_{j}^{\alpha} {\bf S}_{i}^{\beta} ]` is added
for off-diagonal interactions (:math:`\alpha\neq\beta`),


Below is an example of **pair.txt** for 2 site Heisenberg model ::

 0 1 x x 0.5
 0 1 y y 0.5
 0 1 z z 0.5

One can also specify the off-diagonal interaction as ::

 0 1 x x 0.5
 0 1 y y 0.5
 0 1 z z 0.5
 0 1 x y 0.5
 0 1 x z 0.5
 0 1 y z 0.5

Note that interaction :math:`(\alpha,\beta)` and :math:`(\beta,\alpha)` generate the same interaction
since we assume :math:`J_{i,j}^{\alpha,\beta}` is real. 

.. Note that interaction terms must be specified for **(x,y), (x,z), (y,z)** and **(y,x), (z,x), (z,y) cannot be used**.

Exercise
"""""""""""
As a simple exercise, for small system sizes (e.g. 2-site system),  
please try to
confirm that
:math:`H_{xy} = \sum_{i,j} J_{i,j}^{x,y} [{\bf S}_{i}^{x} {\bf S}_{j}^{y}+{\bf S}_{j}^{x} {\bf S}_{i}^{y} ]`,
:math:`H_{yz} = \sum_{i,j} J_{i,j}^{y,z} [{\bf S}_{i}^{y} {\bf S}_{j}^{z}+{\bf S}_{j}^{y} {\bf S}_{i}^{z} ]`, and
:math:`H_{zx} = \sum_{i,j} J_{i,j}^{z,x} [{\bf S}_{i}^{z} {\bf S}_{j}^{x}+{\bf S}_{j}^{z} {\bf S}_{i}^{x} ]`
are equivalent by calculating the ground-state energy and the excited-state energies.

As an advanced exercise,
please try to make input files
for the **Kitaev model** on the honeycomb lattice
by changing **pair.txt** and **input.txt**.
Since the **Kitaev model** can also be treated in Standard mode,
please confirm the results by Expert mode are
consistent with those by Standard mode.

Another example is the **XY model** on the one-dimensional chain.
Although Standard mode only supports the Heisenberg model,
it is easy to treat **XY model** 
by omitting "CoulombInter  coulombinter.def"  and
"Hund  hund.def" in namelist.def generated for the Heisenberg model.
Please confirm the results by Expert mode are
consistent with those by Standard mode.

.. highlight:: none

Lanczos_Step.dat
----------------

| (For the Lanczos method) This file is outputted to show the process
  information for calculating the eigenvector by Lanczos method. An
  example of the file format for the Lanczos method is shown as follows.
| For ``method="Lanczos"``

.. raw:: latex

   \small

::

    LanczosStep  E[1] E[2] E[3] E[4] Target:E[3] E_Max/Nsite
    stp = 2 1.2149586211 14.6560471044 xxxxxxxxxx xxxxxxxxx 0.0000000000 xxxxxxxxx
    stp=4 -5.6626980051 3.1523174817 12.4860778911 21.2322666770 12.4860778911 2.6540333346
    stp=6 -8.5113374325 -2.3219709559 4.3459108959 11.5079386600 4.3459108959 3.0307814358
    stp=8 -9.5061025854 -5.2757708534 -0.1734100333 5.2236216333 -0.1734100333 3.2049774861
    stp=10 -9.7541889139 -6.6054773893 -2.9493235242 1.2364826532 -2.9493235242 3.2686702753
    ...
    stp=84 -10.4543987874 -9.8960493865 -9.7550111859 -9.7407358084 -9.7550111859 3.3731105157
    stp=86 -10.4543987874 -9.8960493865 -9.7550111859 -9.7407358084 -9.7550111859 3.3731105157
    stp=88 -10.4543987874 -9.8960493865 -9.7550111859 -9.7407358084 -9.7550111859 3.3731105157

.. raw:: latex

   \normalsize

| For ``method="CG"``

.. raw:: latex

   \small

::

        Step   Residual-2-norm     Threshold      Energy
            1     6.79819e+00     8.19743e-07    7.86586e+00     8.19743e+00     8.02804e+00
            2     7.47402e+00     3.69905e-07    3.35827e+00     3.63546e+00     3.69905e+00
            3     5.30943e+00     2.44472e-07   -2.44472e+00    -2.23296e+00    -1.95487e+00
            4     4.52737e+00     5.10297e-07   -5.10297e+00    -4.92390e+00    -4.58682e+00
            5     3.66168e+00     7.14105e-07   -7.14105e+00    -6.91226e+00    -6.44532e+00
            6     3.12717e+00     8.27201e-07   -8.27201e+00    -7.93262e+00    -7.44680e+00
      ...
          152     1.05602e-06     1.04544e-06   -1.04544e+01    -9.89605e+00    -9.89605e+00
          153     1.07401e-06     1.04544e-06   -1.04544e+01    -9.89605e+00    -9.89605e+00
          154     9.45018e-07     1.04544e-06   -1.04544e+01    -9.89605e+00    -9.89605e+00

.. raw:: latex

   \normalsize

.. _file_name_5:

File name
~~~~~~~~~

*  ##_Lanczos_Step.dat

## indicates a header defined by [string02] in a ModPara file.

.. _file_format_28:

File format
~~~~~~~~~~~

*  For``method="Lanczos"``

   stp= [int01] [double01]
   [double02] [double03]
   [double04] [double-a]
   [double-b]

*  For ``method="CG"``

   [int01] [double-c]
   [double-d] [double01]
   [double02] ...

.. _parameters_28:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The iteration number of the Lanczos and LOBCG
   method.

*  [double01], [double02],
   [double03], [double04] ...

   **Type :** Double

   **Description :** Eigenvalues computed with the Lanczos or LOBCG
   method (ascending order). Four and ``exct`` eigenvalues are printed
   for the Lanczos and LOBCG method, respectively (in the above case,
   ``exct=3``). While the degenerate eigenstates are printed as a single
   state in the Lanczos method, they are printed separately in LOBCG
   method. In the above case, we can find there is a degeneracy in the
   first excited state.

*  [double-a]

   **Type :** Double

   **Description :** (Only for the Lanczos method) The eigenvalue used
   for the convergence check. It was specified by ``LanczosTarget`` (in
   the above case, ``LanczosTarget=3``).

* [double-b]

   **Type :** Double

   **Description :** (Only for the Lanczos method) The maximum
   eigenvalue divided by the number of sites. It is the lower limit of
   ``LargeValue`` in TPQ method.

*  [double-c]

   **Type :** Double

   **Description :** (Only for LOBCG method) The maximum of the 2-norm
   of each residual vector. It is used for the convergence check.

*  [double-d]

   **Type :** Double

   **Description :** (Only for LOBCG method) The convergence threshold.
   This is obtained by the value specified with ``LanczosEps`` times the
   absolute value of the energy of the ground state.

.. raw:: latex

   \newpage
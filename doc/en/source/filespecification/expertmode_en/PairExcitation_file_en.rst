.. highlight:: none

.. _Subsec:pairexcitation:

PairExcitation file
-------------------

To compute the dynamical correlation function

.. math:: G_n^{O_l,O_r}(z) = \langle \Phi_n | \hat{O}_l (z + E_n - \hat{\cal H})^{-1} \hat{O}_r| \Phi_n \rangle,

we set a pair-excitation operator as

.. math::

    \hat{O}_{l,r} = \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2}
    c_{i \sigma_1}c_{j \sigma_2}^{\dagger} \quad \textrm{or} \quad
    \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2}
    c_{i\sigma_1}^{\dagger}c_{j\sigma_2}

We can compute efficiently by using single :math:`\hat{O}_r` and multiple :math:`\hat{O}_l`.

The type of pair excitation operators (:math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}` or
:math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`) must be the same in the input file.

In the :math:`S_z` conserved system, :math:`\sigma_1` must be equal to :math:`\sigma_2`.

An example of the file format is as follows.

::

    =============================================
    NPair 9
    =============================================
    =============== Pair Excitation =============
    =============================================
    2
    0 0 0 0 1        -0.500000000000000 0.0
    0 1 0 1 1         0.500000000000000 0.0
    2
    0 0 0 0 1        -0.500000000000000 0.0
    0 1 0 1 1         0.500000000000000 0.0
    2
    1 0 1 0 1        -0.500000000000000 0.0
    1 1 1 1 1         0.500000000000000 0.0
    2
    2 0 2 0 1        -0.500000000000000 0.0
    2 1 2 1 1         0.500000000000000 0.0
    2
    3 0 3 0 1        -0.500000000000000 0.0
    3 1 3 1 1         0.500000000000000 0.0
    :

.. _file_format_16:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-:

    Repeat the following block [int01] times.
    The first block corresponds to :math:`\hat{O}_{r}` and other blocks correspond to :math:`\hat{O}_{l}`.
    
    ::
      
        [int02]
        [int03]  [int04]  [int05]  [int06]  [int07]  [double01]  [double02]
        :
        (Repeat [int02] times)

.. _parameters_16:

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of the pair
   excitation operators. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of pair-excitation operators
   :math:`\hat{O}_r` and :math:`\hat{O}_l`.
   For the above example, we have 9 operators (one :math:`\hat{O}_{r}` and 8 :math:`\hat{O}_{l}`),

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of pair-excitation operators
   included in each :math:`\hat{O}_r` and :math:`\hat{O}_l`.

*  [int03], [int05]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int02], [int04] :math:`<` ``Nsite``).

*  [int04], [int06]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | For Hubbard or Kondo system, the index can be selected from
   | 0: Up-spin,
   | 1: Down-spin.
   | For Spin system, the index can be selected from
   | :math:`0, 1, \cdots, 2S+1` (corresponding to
     -:math:`S-0.5, -S+0.5, \cdots S+0.5`).

*  [int06]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a type of pair excitation
     operators:
   | 0: :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger}`
   | 1: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`.

*  [double01], [double02]

   **Type :** Double (a blank parameter is not allowed)

   **Description :** [double01] gives the real part
   of
   :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger} ( c_{i\sigma_1}^{\dagger}c_{j\sigma_2})`,
   while [double02] gives the imaginary part of
   :math:`c_{i\sigma_1}c_{j\sigma_2}^{\dagger} ( c_{i\sigma_1}^{\dagger}c_{j\sigma_2})`.

.. _use_rules_16:

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  A program is terminated when the components of the pair excitation
   operators are double counted.

*  A program is terminated when the conditions :math:`i=j` and
   :math:`k=l` are not satisfied for the spin model.

*  A program is terminated when [int01] is different
   from the total number of two-body Green’s functions defined in this
   file.

*  A program is terminated when
   [int02]-[int06] are outside
   the range of the defined values.

.. raw:: latex

   \newpage

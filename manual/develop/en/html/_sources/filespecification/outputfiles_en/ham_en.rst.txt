.. highlight:: none

.. _Subsec:ham:

ham.dat
-------

| (For the FullDiag method) When ``OutputHam=1`` in the ``CalcMod``
  file, the Hamiltonian calculated by :math:`{\mathcal H}\Phi` is outputted by the
  MatrixMarket format. The recalculation by using this file can be down
  when ``InputHam=1`` in the ``CalcMod`` file. An example of the file
  format is as follows.

::

    %%%%MatrixMarket matrix coordinate complex hermitian
    28 28 56
    1 1 1.000000 0.000000
    2 1 0.500000 0.000000
    3 2 0.500000 0.000000
    4 3 0.500000 0.000000
    5 4 0.500000 0.000000
    6 5 0.500000 0.000000
    7 6 0.500000 0.000000
    7 7 1.000000 0.000000
        …

.. _file_name_15:

File name
~~~~~~~~~

*  ##_ham.dat

## indicates [string02] in a ModPara file.

.. _file_format_39:

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2:
   [int01] [int02] [int03]

*  Lines3-:[int04] [int05] [double01] [double02]

.. _parameters_39:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** Number of rows of Hamiltonian.

*  [int02]

   **Type :** Int

   **Description :** Number of columns of Hamiltonian.

*  [int03]

   **Type :** Int

   **Description :** Number of nonzero elements of Hamiltonian.

*  [double01], [double02]

   **Type :** Double

   **Description :** The value of Hamiltonian;
   [double01] and [double02]
   represent real and imaginary part of the Hamiltonian, respectively.

Sector FullDiag Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``TransSym`` (:ref:`Subsec:TransSym`) and expert-mode ``CalcType=2``,
``OutputHam=1`` writes the :math:`D\times D` Hamiltonian of the selected
sector to ``output/##_ham_sector.mtx`` and stops; no eigenvalues,
correlations or eigenstates are computed. ``InputHam=1`` reads the same file
and diagonalizes that matrix. Such a run can still write, besides the
eigenvalues, the requested correlations and, with ``OutputEigenVec=1``, every
eigenstate. ``OutputHam`` and ``InputHam`` cannot be combined.

The file is a Matrix Market ``coordinate complex hermitian`` matrix. Its basis
is HPhi's ordered sector representatives with their phase convention, and the
``% HPHI key=value`` lines identify that basis. An example for a spin chain
with ``TransSym`` (6 sites, momentum :math:`k=1`, :math:`D=3`):

::

    %%MatrixMarket matrix coordinate complex hermitian
    % HPHI format=sector-hamiltonian-v1
    % HPHI scalar=complex-binary64
    % HPHI model=Spin
    % HPHI nsite=6
    % HPHI nup=3
    % HPHI ndown=3
    % HPHI ne=3
    % HPHI local_site_mask=0000000000000000
    % HPHI fixed_flags=0
    % HPHI phase=1
    % HPHI full_dim=20
    % HPHI sector_dim=3
    % HPHI group_digest=8a08ff3f5c2433e3
    % HPHI sector_count=3
    % HPHI sector_xor=ba8fc1732fc70785
    % HPHI sector_sum=640126bcd43ffa9b
    % HPHI basis_order=110eb94c9a548c63
    3 3 6
    1 1 0.80000000000000004 0
    2 1 0.62 -0.11
    3 1 0.21473720558371184 0.59193575034635193
    2 2 0.48999999999999999 0
    3 2 0.40526279441628799 0.48193575034635183
    3 3 0.48999999999999999 0

*  Line 1 is ``%%MatrixMarket matrix coordinate complex hermitian``.

*  The ``% HPHI`` lines hold each of the following 17 keys once, in any
   order. Integers are decimal; masks and digests have 16 hexadecimal digits:
   ``format=sector-hamiltonian-v1``, ``scalar=complex-binary64``, ``model``
   (``Spin``, ``SpinGC``, ``SpinlessFermion``, ``Hubbard``, ``tJ``,
   ``Kondo``, ``KondoNConserved``, ``KondoGC``, ``HubbardNConserved``,
   ``tJNConserved``), ``nsite``, ``nup``,
   ``ndown``, ``ne``, ``local_site_mask``, ``fixed_flags``, ``phase``,
   ``full_dim``, ``sector_dim``, ``group_digest``, ``sector_count``,
   ``sector_xor``, ``sector_sum`` and ``basis_order``. The values are the
   basis identity of a sector checkpoint; some of them, such as the group and
   sector digests, also appear in the sector manifest of
   :ref:`Subsec:TransSym`.

*  The dimension line is ``D D nnz`` and each data line is
   ``row column real imag``. Only the one-based lower triangle
   (:math:`{\rm row}\ge{\rm column}`) is written. HPhi writes the nonzero
   elements column by column with increasing rows, using 17 significant digits
   (a signed zero is written as +0), so the file restores the binary64 values
   that were diagonalized.

On input the data lines may come in any order and may include explicit zeros.
Missing elements are zero and the upper triangle is the complex conjugate of
the lower one. Diagonal imaginary parts must be exactly zero. Blank lines and
general ``%`` comments may appear between line 1 and the dimension line and
after the last data line; LF and CRLF line ends are accepted.

Identity keys that do not match the sector of the current input (including
another sector of the same dimension), a missing, repeated or unknown key, a
malformed value, a dimension or ``nnz`` mismatch, an index out of range or in
the upper triangle, a repeated coordinate (zeros included), a non-finite
value, a nonzero diagonal imaginary part, fewer or more data lines than
declared, an extra token and an overlong line are errors: every MPI rank
stops before the diagonalization.

Keep the ``% HPHI`` lines when editing coefficients. A matrix made by other
tools is read under the same rules if its basis identity matches; differing
from the Hamiltonian of the definition files is not an error. For example,
this POSIX shell snippet doubles every coefficient and hence every eigenvalue:

::

    awk '/^%/ || NF != 4 { print; next }
         { $3 = sprintf("%.17g", 2 * $3); $4 = sprintf("%.17g", 2 * $4); print }' \
        output/zvo_ham_sector.mtx > edited.mtx
    mv edited.mtx output/zvo_ham_sector.mtx
    # Set OutputHam 0 and InputHam 1 in calc.def and run the same input.

The matrix file of a basis can be read with another MPI rank count and solver
(LAPACK, ScaLAPACK or ELPA). Even with ELPA on several ranks each rank keeps
only its own columns; the root reads the data lines and distributes them in
fixed-size chunks. Eigenstate checkpoints, in contrast, need the rank count
and layout they were written with.

An ``InputHam=1`` run records ``matrix_source=input`` and the ``matrix_digest``
of the values read in the manifest, and its eigenstate checkpoints carry that
digest. TE or CG runs that use these states evolve or re-evaluate them with the
Hamiltonian of their own definition files (:ref:`Subsec:eigenvec`), so their
energies need not equal the eigenvalues of the edited matrix.

.. raw:: latex

   \newpage
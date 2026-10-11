.. highlight:: none

.. _Subsec:eigenvec:

eigenvec.dat
------------

When ``OutputEigenVec=1`` in a ``CalcMod`` file, the eigenvectors are
written in binary form. Iterative solvers write the requested states,
whereas FullDiag writes all eigenvectors.

``InputEigenVec=1`` can read the iterative-solver files for the
corresponding restart or spectrum workflows. FullDiag output uses the
same binary layout, but it is not a FullDiag restart mechanism:
FullDiag always performs the diagonalization.

.. _file_name_18:

File name
~~~~~~~~~

*  ##_eigenvec\_&&\_rank\_$$.dat

## indicates [string02] in a ModPara file, && is the zero-based
eigenstate index, and $$ is a rank number.

For iterative solvers, each MPI rank writes its local Hilbert-space
slice to the file with its own rank number. For FullDiag, each file
contains a complete eigenvector; therefore only
``##_eigenvec_&&_rank_0.dat`` is written, even when the distributed
ScaLAPACK or ELPA solver and ``ExpecMode=1`` or ``2`` are used.

.. _file_format_42:

File format
~~~~~~~~~~~

| This file is written through the following source code (a little
  different fron the actual :math:`{\mathcal H}\Phi` source).

::

    fp = fopen("zvo_eigenvec_0_rank_0.dat", "wb");
    fwrite(&number_of_iterations, sizeof(int), 1,fp);
    fwrite(&local_size, sizeof(unsigned long int),1,fp);
    fwrite(&eigen_vector[0], sizeof(complex double),local_size+1, fp);
    fclose(fp);

where ``number_of_iterations`` is the number of iterations,
``local_size`` is the size of eigenvector (if MPI is used, it differs
from the dimension of the Hilbert space), ``eigen_vector`` is the
(complex) eigenvector.

**Note:** The first component of ``eigen_vector`` (``eigen_vector[0]``)
is not used for calculation.

For FullDiag, ``number_of_iterations`` is zero and ``local_size`` is
the full Hilbert-space dimension. Thus, enabling this option writes
:math:`N_{\rm H}` files containing :math:`N_{\rm H}` complex
components each; disk usage is :math:`O(N_{\rm H}^2)`.

Sector FullDiag eigenstates
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``TransSym`` in expert mode, ``CalcType=2`` and ``OutputEigenVec=1``
save **every** sector eigenstate as
``##_eigenvec_&&_rank_$$.dat``. The state number && is zero-based in
ascending energy, the same order as ``##_energy_sector.dat``. States within
a degenerate level are not unique: their order, phase and basis can differ
between solvers or runs. Every MPI rank writes its own block, including
ranks that own no component, so a sector of dimension :math:`D` run on
:math:`P` ranks writes :math:`D\times P` files (160,000 files for
:math:`D=10{,}000` and :math:`P=16`). Selecting a range of states is not
available.

These files use the sector checkpoint format described in
:ref:`Subsec:TransSym`, not the raw binary layout above. The names are
deliberately the same as for raw FullDiag; whether ``TransSym`` is used
decides the format, so check the run before reading a file. The header
records ``source_method=2`` (FullDiag), the state number, step 0, time 0 and
the ``hphi-sector-matrix-fnv1a64-v1`` digest of the diagonalized matrix
instead of the definition-file Hamiltonian digest.

A saved state can be reused in the same sector, with the same MPI rank
count, basis order, phase and layout:

* Time evolution: copy the rank files of the chosen state, for example
  state 2, to ``output/seed_rank_$$.dat`` of a separate run directory and
  add ``SpectrumVec seed`` to the namelist (``CalcType=4``,
  ``InputEigenVec=1``). TE writes its own ``##_eigenvec_&&`` files, which
  would overwrite other FullDiag states in the same directory.
* CG: ``InputEigenVec=1`` reads ``##_eigenvec_0`` to
  ``##_eigenvec_(exct-1)`` and evaluates them with the CG input Hamiltonian
  without iterating. Use ``OutputEigenVec=0`` for this re-evaluation;
  ``OutputEigenVec=1`` overwrites the input files with CG checkpoints. Copy
  the FullDiag files to a separate directory to keep the originals.

TE and CG use ``HPHI_SYMMETRY_BASIS_LAYOUT=distributed`` by default, while
FullDiag writes replicated-layout states. Set
``HPHI_SYMMETRY_BASIS_LAYOUT=replicated`` for the consumer run. The consumer
evaluates the Hamiltonian of its own definition files; a quench is allowed.
Because the FullDiag digest is a matrix digest, the load log reports
``source_hamiltonian_kind=sector_matrix_v1 hamiltonian_changed=unknown``.

A portable POSIX shell example for a TE run from state 2, after a two-rank
FullDiag run in ``fulldiag/``:

::

    mkdir -p te/output
    cp fulldiag/*.def te/
    for rank in 0 1; do
      cp fulldiag/output/zvo_eigenvec_2_rank_${rank}.dat \
         te/output/seed_rank_${rank}.dat
    done
    # In te/: set CalcType 4 and InputEigenVec 1 in calc.def, add a
    # TEOneBody time grid and "SpectrumVec seed" to the namelist, then run
    # with the same two MPI ranks:
    # HPHI_SYMMETRY_BASIS_LAYOUT=replicated mpiexec -np 2 HPhi -e namelist.def

.. raw:: latex

   \newpage

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

.. raw:: latex

   \newpage

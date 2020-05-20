.. highlight:: none

.. _Subsec:eigenvec:

eigenvec.dat
------------

When OutputEigenVec=1 in a CalcMod file, the eigenvectors calculated by
the Lanczos method are outputted. When InputEigenVec=1 in a CalcMod
file, eigenvectors are inputted by this outputted file. The file format
is of the binary type.

.. _file_name_18:

File name
~~~~~~~~~

*  ##_eigenvec\_&&\_rank\_$$.dat

## indicates [string02] in a ModPara file, && is the number of
eigenvalues, and $$ is a number of rank.

.. _file_format_42:

File format
~~~~~~~~~~~

| This file is written through the following source code (a little
  different fron the actual :math:`{\mathcal H}\Phi` source).

::

    fp = fopen("zvo_eigenvec_0_rank_0.dat", "wb");
    fwrite(&number_of_interations, sizeof(int), 1,fp);
    fwrite(&local_size, sizeof(unsigned long int),1,fp);
    fwrite(&eigen_vector[0], sizeof(complex double),local_size+1, fp);
    fclose(fp);

where ``number_of_interations`` is the number of iterations,
``local_size`` is the size of eigenvector (if MPI is used, it differs
from the dimension of the Hilbert space), ``eigen_vector`` is the
(complex) eigenvector.

**Note:** The fist component of ``eigen_vector`` (``eigen_vector[0]``)
is not used for calculation.

.. raw:: latex

   \newpage
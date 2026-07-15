.. highlight:: none

.. _Subsec:restart:

tmpvec.dat
----------

When ReStart=1, 2 in a CalcMod file, vectors after the calculation stops
at an indicated step are outputted. The file format is of the binary
type. An example of the file format is as follows.

.. _file_name_19:

File name
~~~~~~~~~

*  Lanczos method: ##_tmpvec_rank_$$.dat

*  TPQ and LOBPCG method: ##_tmpvec_set_&&_rank_$$.dat .

## indicates [string02] in a ModPara file, and $$ is the number of rank.
&& is the sampling number for the TPQ calculation.

.. _file_format_43:

File format
~~~~~~~~~~~

This file is written through the following source code (a little
different from the actual :math:`{\mathcal H}\Phi` source).

| Lanczos

::

    fp = fopen("zvo_tmpvec_rank_0.dat", "wb");
    fwrite(&number_of_interations, sizeof(int), 1,fp);
    fwrite(&local_size, sizeof(unsigned long int),1,fp);
    fwrite(&last_vector[0], sizeof(complex double),local_size+1, fp);
    fwrite(&second_last_vector[0], sizeof(complex double),local_size+1, fp);
    fclose(fp);

| TPQ and LOBPCG

::

    fp = fopen("zvo_tmpvec_set_0_rank_0.dat", "wb");
    fwrite(&number_of_interations, sizeof(int), 1,fp);
    fwrite(&local_size, sizeof(unsigned long int),1,fp);
    fwrite(&last_vector[0], sizeof(complex double),local_size+1, fp);
    fclose(fp);

where ``number_of_interations`` is the number of iterations,
``local_size`` is the size of eigenvector (if MPI is used, it differs
from the dimension of the Hilbert space), ``last_vector`` is the vector
at the last iteration and ``second_last_vector`` is the vector at the
second last iteration.

**Note:** The fist component of ``last_vector`` and
``second_last_vector`` (``last_vector[0]`` and
``second_last_vector[0]``) are not used for calculation.

.. raw:: latex

   \newpage

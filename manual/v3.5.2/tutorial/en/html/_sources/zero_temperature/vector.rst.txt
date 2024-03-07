Use eigenvectors
^^^^^^^^^^^^^^^^^^^^^^^^^
In this tutorial, we will study how to read the eigenvectors.
In the standard mode, setting ``EigenVecIO = "Out"`` makes HPhi write the calculated eigenvectors as ``output/zvo_eigenvec_[index]_rank_[rank].dat``, where ``[index]`` is the index of the states (e.g., the ground state has ``[index] = 0``) and ``[rank]`` is the rank of the process.
In the MPI parallelization with :math:`N_{\text{para}}` processes, HPhi splits the whole Hilbert space into the :math:`N_{\text{para}}` blocks and each process treats one of them.
The file format is described in the `reference manual <http://issp-center-dev.github.io/HPhi/manual/master/en/html/filespecification/outputfiles_en/tmpvec_en.html>`_ .
For example, the following python function (``samples/tutorial-1.7/read_gs.py``) reads the vector::

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
          ret = np.zeros(nelems, dtype=np.complex128)
          f.read(16)
          for i in range(nelems):
              re = unpack("d", f.read(8))[0]
              im = unpack("d", f.read(8))[0]
              ret[i] = complex(re, im)
          return ret

Exercise
"""""""""""
Check the orthogonality of the eigenvectors calculated by the LOBCG method by calculating the norm and the inner product of some of the eigenvectors.

*Hint* : In the standard mode, the ``exct`` keyword controls the number of eigenvectors to be calculated.

*Solution* : See ``samples/tutorial-1.7/solution.py``.
This script firstly generates the input file for calculating the ground state and the first excited state of the :math:`L=8` AFH chain with the LOBCG method, next invokes ``HPhi``, then reads the vectors, and finally calculates the norms and the inner product.

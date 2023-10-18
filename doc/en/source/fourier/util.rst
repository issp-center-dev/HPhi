Behavior of ``greenr2k`` utility
================================

This utility is used as follows:

.. code-block:: bash

   $ ${PATH}/greenr2k ${NAMELIST} ${GEOMETRY}

where ``${PATH}`` is the path to the directory where
the executable ``greenr2k`` exists,
${NAMELIST} is the NameList input-file name of :math:`{\mathcal H}\Phi`/mVMC, and
${GEOMETRY} is the path to the :ref:`geometry` file.

The behavior of this utility is slightly different between the correlation functions from
each mode of :math:`{\mathcal H}\Phi` (Lanczos, TPQ, Full diagonalization, LOBCG)
and mVMC.
In the following cases, we assume that
``CDataFileHead`` in the ModPara input file is ``"zvo"`` (default).

HPhi-Lanczos
~~~~~~~~~~~~

In this case, ``HPhi`` writes correlation functions to the files
``zvo_cisajs.dat`` (one body) and ``zvo_cisajscktalt.dat`` (two body)
in ``output/`` directory.
``greenr2k`` utility reads this files, performs the Fourier transformation, and
generate single file ``zvo_corr.dat`` in ``output/`` directory.

HPhi-TPQ
~~~~~~~~

``HPhi`` writes correlation functions to files
``zvo_cisajs_run*step*.dat`` (one body), ``zvo_cisajscktalt_run*step*.dat`` (two body)
at each trial and TPQ step to the ``output/`` directory.
``greenr2k`` utility reads the one- and the two-body correlation function at each trial/TPQ-step,
and performs Fourier transformation, and
write to a file ``zvo_corr_run*step*.dat`` in ``output/`` directory.

HPhi-Full diagonalization and LOBCG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``HPhi`` writes correlation functions to files
``zvo_cisajs_eigen*.dat`` (one body) and ``zvo_cisajscktalt_eigen*.dat`` (two body)
for each wavefunction to the ``output/`` directory.
``greenr2k`` utility reads the one- and the two-body correlation function at each state
and performs Fourier transformation, and
write to a file ``zvo_corr_eigen*.dat`` in ``output/``.

If the temperature which is the optional imput is specified,
the temperature-dependent correlation function is written into files ``zvo_corr_t*.dat`` in the ``output/`` directory.

mVMC
~~~~

``vmc.out`` performs calculations according to the input parameters ``NDataIdxStart`` and ``NDataQtySmp``
in ``ModPara`` file, and it generates
``zvo_cisajs_???.dat`` (one body) and ``zvo_cisajscktalt_???.dat`` (two body)
in ``output/`` directory.
``greenr2k`` utility reads all of these files, performs Fourier transformation,
computes the average 

.. math::

   \begin{align}
   \langle A \rangle = \frac{1}{N_{\rm Try}} \sum_{i=1}^{N_{\rm Try}} A_i
   \end{align}

and the standard error

.. math::
   
   \begin{align}
   \delta A = \frac{1}{N_{\rm Try} - 1}
   \sqrt{\frac{1}{N_{\rm Try}} \sum_{i=1}^{N_{\rm Try}} (A_i - \langle A \rangle)^2}
   \end{align}

of the real- and imaginary-part of each correlation function, and
writes them to a file ``zvo_corr_eigen*.dat`` in ``output/`` directory.

Behavior of ``dynamicalr2k`` utility
====================================

This utility is used as follows:

.. code-block:: bash

   $ ${PATH}/dynamicalr2k ${NAMELIST} ${GEOMETRY}

where ``${PATH}`` is the path to the directory where
the executable ``dynamicalr2k`` exists,
${NAMELIST} is the NameList input-file name of :math:`{\mathcal H}\Phi`, and
${GEOMETRY} is the path to the :ref:`geometry` file.

In the default setting, ``HPhi`` generate site-dependent dynamical green's function into files ``zvo_DynamicalGreen_*.dat`` in ``output/`` directory.
``dynamicalr2k`` reads these files, performes Fourier transform, and output files ``zvo_dyn*.dat`` in ``output/`` directory.

If the temperature which is the optional parameter is specified, the temperature-dependent Green's funtion is outputted into files ``zvo_dyn_t*.dat`` in ``output/`` directory.

.. highlight:: none

Full Diagonalization method
===========================

Overview
--------

We generate the matrix of :math:`\hat{\mathcal H }` by using the real space configuration 
:math:`| \psi_j \rangle`\(:math:`j=1\cdots d_{\rm H}`, where :math:`d_{\rm H}` is the dimension of the Hilbert space): 
:math:`{\mathcal H }_{ij}= \langle \psi_i | \hat {\mathcal H } | \psi_j \rangle`.
By diagonalizing this matrix,
we can obtain all the eigenvalues :math:`E_{i}` and eigenvectors :math:`|\Phi_i\rangle` (:math:`i=1 \cdots d_{\rm H}`). 
In the diagonalization, we use a LAPACK routine, such as ``dsyev`` or ``zheev``.
We also calculate and output
the expectation values :math:`A_i \equiv \langle \Phi_i | {\hat A} | \Phi_i\rangle`.
These values are used for the finite-temperature calculations.

Finite-temperature calculations
-------------------------------

From
:math:`A_i \equiv \langle \Phi_i | {\hat A} | \Phi_i\rangle`,
we calculate the finite-temperature properties by using the relation 

.. math::

   \langle {\hat A}\rangle=\frac{\sum_{i=1}^N A_i {\rm  e}^{-\beta E_i}}{\sum_{i=1}^N{\rm  e}^{-\beta E_i}}.

The calculation should be performed by using the own postscripts.

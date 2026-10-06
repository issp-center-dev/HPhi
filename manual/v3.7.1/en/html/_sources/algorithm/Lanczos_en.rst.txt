.. highlight:: none

Lanczos method
==============

Details of Lanczos method
-------------------------

Some parts of this section are based on the manual of TITPACK [#]_ and the textbook published by M. Sugihara and K. Murota [#]_ (these references are written in Japanese).

In the Lanczos method, by successively operating the Hamiltonian 
to the initial vector, we obtain the accurate eigenvalues around
the maximum and minimum eigenvalues and associated eigenvectors.  
Because we can perform the Lanczos method by using only two
vectors, the dimensions of which are the dimensions of the total Hilbert space [#]_ , the Lanczos method is frequently used for the 
diagonalization of the large matrices.
As explained in detail below,
one additional vector is necessary for
obtaining the eigenvector.

The principle of the Lanczos method is
based on the power method.
In the power method,
by successively operating the Hamiltonian :math:`\hat{\mathcal H }` to the
arbitrary vector :math:`\boldsymbol{x}_{0}`, we generate :math:`\hat{\mathcal H }^{n}\boldsymbol{x}_{0}`.
The obtained space 
:math:`\mathcal{K}_{n+1}(\hat{\mathcal H },\boldsymbol{x}_{0})=\{\boldsymbol{x}_{0},\hat{\mathcal H }^{1}\boldsymbol{x}_{0},\dots,\hat{\mathcal H }^{n}\boldsymbol{x}_{0}\}`
is called the Krylov subspace.
The initial vector is represented by the superposition 
of the eigenvectors 
:math:`\boldsymbol{e}_{i}` (the corresponding eigenvalues are :math:`E_{i}`) of :math:`\hat{\mathcal H }` as 

.. math::

   \boldsymbol{x}_{0}=\sum_{i}a_{i}\boldsymbol{e}_{i}.
   
Here, :math:`E_{0}` denotes the maximum absolute values of the eigenvalues.
We note that all the eigenvalues are real numbers because the Hamiltonian is Hermitian.
By operating :math:`\hat{\mathcal H }^{n}` to the initial vector,
we obtain the relation as

.. math::

   \hat{\mathcal H }^{n}\boldsymbol{x}_{0}=E_{0}^{n}\Big[ a_{0}\boldsymbol{e}_{0}+\sum_{i\neq0}\left(\frac{E_{i}}{E_{0}}\right)^na_{i}\boldsymbol{e}_{i}\Big].

This relation indicates that
the eigenvector of :math:`E_{0}` becomes dominant for sufficiently large :math:`n`. 
In the Lanczos method,
we obtain the eigenvalues and eigenvectors 
by performing the appropriate transformation for the obtained Krylov subspace.

In the Lanczos method,
we successively generate the normalized orthogonal basis 
:math:`\boldsymbol{v}_{0},\dots,\boldsymbol{v}_{n-1}` from the Krylov subspace :math:`\mathcal{K}_{n}(\hat{\mathcal H },\boldsymbol{x}_{0})`.
We define an initial vector and associated components as 
:math:`\boldsymbol{v}_{0} =\boldsymbol{x}_{0}/|\boldsymbol{x}_{0}|`,
:math:`\beta_{0}=0,\boldsymbol{x}_{-1}=0`.
From this initial condition,
we can obtain the normalized orthogonal basis:

.. math::

   \alpha_{k} &= (\hat{\mathcal H }\boldsymbol{v}_{k},\boldsymbol{v}_{k}), \\
   \boldsymbol{w}   &= \hat{\mathcal H }\boldsymbol{v}_{k}-\beta_{k}\boldsymbol{v}_{k-1}-\alpha_{k}\boldsymbol{v}_{k}, \\
   \beta_{k+1} &= |\boldsymbol{w}|, \\
   \boldsymbol{v}_{k+1} &= \frac{\boldsymbol{v}_{k}}{|\boldsymbol{v}_{k}|}.\\

From these definitions, it it obvious that :math:`\alpha_{k}`, :math:`\beta_{k}` are real numbers.

In the subspace spanned by these normalized orthogonal basis,
the Hamiltonian is transformed as

.. math::

   T_{n}=V_{n}^{\dagger}\hat{\mathcal H } V_{n}.

Here,
:math:`V_{n}` is a matrix whose column vectors are :math:`\boldsymbol{v}_{i}(i=0,1,\dots,n-1)`.
:math:`T_{n}` is a tridiagonal matrix and its diagonal elements
are :math:`\alpha_{i}` and
subdiagonal elements are :math:`\beta_{i}`.
It is known that
the eigenvalues of :math:`\hat{\mathcal H }` are well approximated by 
the eigenvalues of :math:`T_{n}` for sufficiently large :math:`n`.
(We note that :math:`V^{\dagger}V=I`, :math:`I` is an identity matrix).
The original eigenvectors of :math:`\hat{\mathcal H }` are obtained 
by :math:`\boldsymbol{e}_{i}=V\tilde{\boldsymbol{e}}_{i}`,
where  :math:`\tilde{\boldsymbol{e}}_{i}` denotes
the eigenvectors of :math:`T_{n}`.
From :math:`V`, 
we can obtain the eigenvectors of :math:`\hat{\mathcal H }`
by performing the Lanczos method.
However, in the actual calculations,
it is difficult to keep :math:`V`, because its dimension
is large [dimension of :math:`V` = (dimension of the total Hilbert space) :math:`\times` (the number of Lanczos iterations)].
Thus, to obtain the eigenvectors, 
we again perform the same Lanczos calculations
after we obtain the eigenvalues from the Lanczos methods.
In the first Lanczos calculation, we keep :math:`\tilde{\boldsymbol{e}_{i}}`, 
because its dimension is small [#]_ .
From this procedure, we obtain the eigenvectors  from :math:`V`.

In the Lanczos method,
within a few hundred or thousand Lanczos iterations,
we obtain accurate eigenvalues near the maximum and minimum eigenvalues.
The necessary number of iterations is sufficiently small as 
compared to the dimensions
of the total Hilbert space.

We note that it is shown that
the errors of the maximum and minimum eigenvalues
become exponentially small as a function of Lanczos iterations 
(for details, see Ref. [2]_ ).

Inverse iteration method
------------------------

From the approximate value of the eigenvalues :math:`(E_{n})`,
by successively operating :math:`(\hat{\mathcal H }-E_{n})^{-1}`
to the initial vector :math:`\boldsymbol{y}_{0}`,
we can obtain the accurate eigenvector for :math:`E_{n}`.

From :math:`(\hat{\mathcal H }-E_{n})^{-1}\boldsymbol{y}_{0}`,
we obtain linear simultaneous equations such as  

.. math::

   \boldsymbol{y}_{k}=(\hat{\mathcal H }-E_{n})\boldsymbol{y}_{k+1}.

By solving this equation using the
conjugate gradient method (CG method),
we obtain the eigenvector.
From the obtained eigenvector,
we can calculate the eigenvalues and correlation functions. 
We note that four additional vectors are necessary to
perform the CG method.
For a large system size,
it may be impossible to allocate memory to the
additional vectors.

Details of implementation
-------------------------

**Initial vector**
^^^^^^^^^^^^^^^^^^

For the Lanczos method, an initial vector is specified with ``initial_iv``:math:`(\equiv r_s)` defined in an input file for Standard mode or a ModPara file for Expert mode. The type of initial vector can be selected as a real number or complex number by using ``InitialVecType`` in a ModPara file.


 * For canonical ensemble and ``initial_iv``:math:`\geq 0`,
   a component of a target of the Hilbert space is given by
   
   .. math::
 
     (N_{\rm dim}/2 + r_s ) \% N_{\rm dim},

   where :math:`N_{\rm dim}` is the total number of the Hilbert spaces and :math:`N_{\rm dim}/2` is added to avoid selecting a special Hilbert space for a default value ``initial_iv`` :math:`=1`.
   When the type of initial vector is selected as a real number, the coefficient value is given by :math:`1`, while when it is selected as a complex number, the value is given by :math:`(1+i)/\sqrt{2}`.

 * For a grand canonical ensemble or ``initial_iv`` :math:`<0`,
   the initial vector is given by using a random generator, i.e., the coefficients of all the components for the initial vector are given by random numbers. The seed is calculated as 
   
   .. math::
   
      123432+|r_s|,

   where :math:`r_s` is the number given by an input file and :math:`n_{\rm run}` is the number of runs. The maximum value of :math:`n_{\rm run}` is defined by ``NumAve`` in an input file for Standard mode or a ModPara file for Expert mode. Random numbers are generated by using SIMD-oriented Fast Mersenne Twister (dSFMT) [#]_ . 

**Convergence condition**
^^^^^^^^^^^^^^^^^^^^^^^^^

In :math:`{\mathcal H}\Phi`,
we use ``dsyev`` (routine of LAPACK)
for diagonalization of :math:`T_{n}`.
We use the energy of the first excited state of :math:`T_{n}`
as the criterion of convergence. 
In the standard setting,
after five Lanczos steps,
we diagonalize :math:`T_{n}` every two Lanczos steps.
If the energy of the first excited states coincides with
the previous energy within the specified accuracy,
the Lanczos iteration finishes.
The accuracy of the convergence can be specified by 
``CDataFileHead``\(ModPara file in the expert mode).

After obtaining the eigenvalues,
we again perform the Lanczos iteration
to obtain the eigenvector.
From the eigenvectors :math:`|n\rangle`,
we calculate 
energy :math:`E_{n}=\langle n|\hat{\mathcal H }|n\rangle` and
variance :math:`\Delta=\langle n|\hat{\mathcal H }^{2}|n\rangle -(\langle n|\hat{\mathcal H }|n\rangle)^2`.
If :math:`E_{n}` coincides with the eigenvalues obtained by the Lanczos iteration and 
:math:`\Delta` is smaller than the specified value,
we finish diagonalization.

If the accuracy of the Lanczos method is not sufficient,
we perform the CG method to obtain the eigenvector.
As an initial vector of the CG method,
we use the eigenvectors obtained by the Lanczos method in the standard setting.
This frequently accelerates the convergence.

.. [#] \http://www.stat.phys.titech.ac.jp/~nishimori/titpack2_new/index-e.html
.. [#] \M. Sugihara, K. Murota, Theoretical Numerical Linear Algebra, Iwanami Stud-ies in Advanced Mathematics, Iwanami Shoten, Publishers, 2009.
.. [#] \In :math:`{\mathcal H}\Phi`, to reduce the numerical cost, we use some additional vectors; a vector for accumulating the real-space diagonal elements of the Hamiltonian and a vector for specifying the given :math:`S_{z}` space and given particle space. The dimension of these vectors is that of the Hilbert space.
.. [#] \Upper bound of the dimensions of :math:`\tilde{\boldsymbol{e}_{i}}` is \# of Lanczos iterations.
.. [#] \http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html

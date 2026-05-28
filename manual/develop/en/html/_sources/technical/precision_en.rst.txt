Numerical Precision
===================

This section describes numerical precision considerations in :math:`{\mathcal H}\Phi`.

Floating-Point Representation
-----------------------------

:math:`{\mathcal H}\Phi` uses IEEE 754 double precision floating-point numbers throughout:

- Real numbers: 64-bit double (``double``)
- Complex numbers: 128-bit complex double (``double complex``)

This provides approximately 15-16 significant decimal digits of precision.

Sources of Numerical Error
--------------------------

**Finite precision arithmetic**
  Each floating-point operation introduces small rounding errors.
  These accumulate during iterative algorithms.

**Cancellation errors**
  Subtracting nearly equal numbers can lose significant digits.
  This can occur in certain physical quantities near phase transitions.

**Overflow and underflow**
  Very large or small intermediate values can exceed the representable range.
  :math:`{\mathcal H}\Phi` is designed to avoid these in normal usage.

Convergence Criteria
--------------------

:math:`{\mathcal H}\Phi` uses several convergence criteria:

**Lanczos eigenvalue convergence**
  The Lanczos iteration converges when the relative change in eigenvalue
  falls below a specified threshold (default: :math:`10^{-14}`, i.e.
  ``LanczosEps = 14``).

**Eigenvector residual**
  The eigenvector accuracy is checked via :math:`|H|\psi\rangle - E|\psi\rangle|`.

**Orthogonality**
  For multiple eigenvalues, orthogonality between eigenvectors is maintained
  through reorthogonalization.

Setting Convergence Thresholds
------------------------------

The convergence thresholds are specified in the ``ModPara`` input file:

``LanczosEps``
  Integer controlling the convergence threshold. For the Lanczos method,
  the iteration is judged to have converged when the relative error
  between an eigenvalue and that of the previous step becomes smaller
  than :math:`10^{- {\tt LanczosEps}}` (default: ``LanczosEps = 14``).
  For ``method="CG"``, the calculation finishes when the 2-norm of the
  residual vector becomes smaller than :math:`10^{- {\tt LanczosEps}/2}`.

``LanczosTarget``
  Integer giving the target eigenvalue used to judge convergence
  (e.g. ``1`` for the ground state, ``2`` for the first excited state).

``Lanczos_max``
  Maximum number of iteration steps; the calculation finishes earlier
  once the threshold above is reached.

For most applications, the default values provide sufficient accuracy.

Validation and Testing
----------------------

:math:`{\mathcal H}\Phi` includes test cases that verify:

- Eigenvalue accuracy against known analytical results
- Hermiticity of the Hamiltonian
- Conservation of particle number and :math:`S_z`
- Consistency between different calculation methods

Recommendations
---------------

**Check physical quantities**
  Always verify that conserved quantities (particle number, :math:`S_z`, etc.)
  remain constant to within numerical precision.

**Compare methods**
  When possible, compare results from Lanczos and LOBPCG methods.
  Agreement indicates reliable results.

**System size dependence**
  Be cautious about extrapolating results to larger systems,
  as numerical errors may scale differently.

**Report anomalies**
  If you observe unexpected numerical behavior, please report it
  to the developers with a minimal reproducing example.

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
  The Lanczos iteration converges when the change in eigenvalue
  falls below a specified threshold (default: :math:`10^{-12}`).

**Eigenvector residual**
  The eigenvector accuracy is checked via :math:`|H|\psi\rangle - E|\psi\rangle|`.

**Orthogonality**
  For multiple eigenvalues, orthogonality between eigenvectors is maintained
  through reorthogonalization.

Setting Convergence Thresholds
------------------------------

In the ``CalcMod`` input file:

``CDataFileHead``
  Prefix for output files

The convergence threshold is typically set in the source code.
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

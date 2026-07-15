.. highlight:: none

.. _Sec:troubleshooting:

Troubleshooting
===============

This section explains error messages that may occur during :math:`{\mathcal H}\Phi` execution and how to resolve them.

Quick Reference
---------------

**Command-line and File Errors**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - Error Message
     - Cause
     - Solution
   * - ``argument count is incorrect``
     - Invalid number of command-line arguments
     - Run as ``HPhi -e namelist.def`` or ``HPhi -s stan.in``
   * - ``Error: The file (filename) does not exist.``
     - Specified file does not exist
     - Check file path and filename
   * - ``ERROR: Open input file``
     - Cannot open input file
     - Verify file existence and read permissions

**Keyword and Parameter Errors**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - Error Message
     - Cause
     - Solution
   * - ``ERROR ! Unsupported Keyword !``
     - Non-existent keyword specified
     - Check manual for supported keywords
   * - ``ERROR ! Keyword (name) is duplicated !``
     - Same keyword specified twice
     - Remove duplicate keyword
   * - ``ERROR ! (name) is NOT specified !``
     - Required keyword not specified
     - Add the required keyword

**Model and Solver Errors**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - Error Message
     - Cause
     - Solution
   * - ``ERROR ! Unsupported Solver : (name)``
     - Unsupported solver specified
     - Use Lanczos, TPQ, FullDiag, CG, or TimeEvolution
   * - ``ERROR ! Unsupported Model : (name)``
     - Unsupported model specified
     - Use Hubbard, Spin, Kondo, or SpinlessFermion
   * - ``Sorry, this system is unsupported in the STANDARD MODE...``
     - Unsupported model/lattice combination
     - Use Expert mode or choose a supported combination

**Memory Errors**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - Error Message
     - Cause
     - Solution
   * - ``Malloc ERROR for XXX``
     - Memory allocation failed
     - Reduce system size or use a machine with more memory
   * - ``There is not enough memory. Please reduce MPI process number.``
     - Not enough memory per process
     - Reduce MPI process count to increase memory per process

**MPI Errors**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30
   :align: left

   * - Error Message
     - Cause
     - Solution
   * - ``You can not use OpenMP and MPI at the same time.``
     - OpenMP and MPI combination is restricted
     - Use either OpenMP or MPI
   * - ``ERROR: MPI process count is not a power of 2``
     - MPI process count is not a power of 2
     - Set MPI process count to 1, 2, 4, 8, 16, etc.
   * - ``ERROR: Not enough MPI processes.``
     - Insufficient MPI processes
     - Increase the process count

**Hamiltonian Validation Errors (Canonical Ensemble)**

The following errors occur when the combination of site count, electron count, and spin is physically impossible in canonical ensemble calculations.

.. list-table::
   :header-rows: 1
   :widths: 50 50
   :align: left

   * - Error Message
     - Meaning
   * - ``ERROR ! abs(2 * Sz) > nsite in Hubbard model !``
     - Absolute value of Sz exceeds half the number of sites
   * - ``ERROR ! Nelec > 2 * nsite in Hubbard model !``
     - Electron count exceeds twice the number of sites
   * - ``ERROR ! (nelec + 2 * Sz) % 2 != 0 in Hubbard model !``
     - Parity of electron count and Sz does not match
   * - ``ERROR ! nelec <= nsite && 2 * |Sz| > nelec in Hubbard model !``
     - Sz is too large for the electron count
   * - ``ERROR ! nelec > nsite && 2 * |Sz| > 2 * nsite - nelec in Hubbard model !``
     - Sz is too large when electron count is high
   * - ``ERROR ! abs(2 * Sz) > nsite in Spin model !``
     - Sz exceeds half the number of sites in spin model
   * - ``ERROR ! (nsite + 2 * Sz) % 2 != 0 in Spin model !``
     - Parity of site count and Sz does not match in spin model
   * - ``ERROR ! Nelec_cond / 2 + Nelec_loc > nsite in Kondo model !``
     - Electron count exceeds site count in Kondo model
   * - ``ERROR ! (nelec_cond + nelec_loc + 2 * Sz) % 2 != 0 in Kondo model !``
     - Parity of electron count and Sz does not match in Kondo model

**Solution:** Check input parameters ``Nsite``, ``Nelec``, ``2Sz`` and correct to physically valid combinations.

**Warning Messages**

.. list-table::
   :header-rows: 1
   :widths: 50 50
   :align: left

   * - Message
     - Description
   * - ``Check ! (name) is SPECIFIED but will NOT be USED.``
     - Specified parameter is not used in current calculation. Comment out if unnecessary
   * - ``(name) = (value) ###### DEFAULT VALUE IS USED ######``
     - Default value is being used (not an error)

----

Detailed Solutions for Common Errors
------------------------------------

Out of Memory Error
^^^^^^^^^^^^^^^^^^^

**Symptom:** ``Malloc ERROR`` is displayed and calculation cannot start.

**Cause:** Hilbert space dimension is too large and required memory cannot be allocated.

**Solutions:**

1. **Reduce system size**

   Reducing the number of sites exponentially decreases the Hilbert space dimension.

2. **Increase MPI process count**

   Using more processes like ``mpirun -np 16 HPhi ...`` distributes the Hilbert space across processes.

3. **Check available memory**

   Review ``CHECK_Memory.dat`` to understand required memory.

MPI Process Count Error
^^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** ``ERROR: MPI process count is not a power of 2``

**Cause:** :math:`{\mathcal H}\Phi` requires MPI process count to be a power of 2 (1, 2, 4, 8, 16, 32, ...).

**Solution:**

.. code-block:: bash

   # Correct examples
   mpirun -np 4 HPhi -e namelist.def
   mpirun -np 8 HPhi -e namelist.def
   mpirun -np 16 HPhi -e namelist.def

   # Wrong example (3 is not a power of 2)
   mpirun -np 3 HPhi -e namelist.def  # Error

File Not Found Error
^^^^^^^^^^^^^^^^^^^^

**Symptom:** ``Error: The file (filename) does not exist.``

**Solutions:**

1. Verify filenames in namelist.def
2. For relative paths, check position relative to execution directory
3. Check case sensitivity of filenames (Linux/macOS are case-sensitive)

Convergence Issues
^^^^^^^^^^^^^^^^^^

**Symptom:** Lanczos or CG method does not converge.

**Solutions:**

1. **Increase Lanczos steps**

   Adjust ``LanczosEps`` or ``Lanczos_max`` in ModPara file

2. **Change initial vector**

   Try different random seeds with ``initial_iv`` parameter

3. **Degenerate eigenvalues**

   Convergence may be slow when degenerate eigenvalues exist

Error messages
==============

*  ``ERROR ! Unsupported Keyword !``

   The program stops because unsupported keyword is specified.

*  ``"ERROR !  Keyword`` *keyword* ``is duplicated !``

   The program stops because a parameter is specified twice.

*  ``ERROR ! Unsupported Solver :`` *solver*

*  ``ERROR ! Unsupported Model :`` *model*

*  | ``Sorry, this system is unsupported in the STANDARD MODE...``
   | ``Please use the EXPART MODE, or write a NEW FUNCTION and post it us.``

   The program stops because unsupported parameter for ``method``,
   ``model``, or ``lattice`` is specified.

*  ``ERROR ! abs(2 * Sz) > nsite in Hubbard model !``

*  ``ERROR ! Nelec > 2 * nsite in Hubbard model !``

*  ``ERROR ! (nelec + 2 * Sz) % 2 != 0 in Hubbard model !``

*  ``ERROR ! nelec <= nsite && 2 * |Sz| > nelec in Hubbard model !``

*  ``ERROR ! nelec > nsite && 2 * |Sz| > 2 * nsite - nelec in Hubbard model !``

*  ``ERROR ! abs(2 * Sz) > nsite in Spin model !``

*  ``ERROR ! (nsite + 2 * Sz) % 2 != 0 in Spin model !``

*  ``ERROR ! abs(2 * Sz) > nsite in Hubbard model !``

*  ``ERROR ! Nelec_cond / 2 + Nelec_loc > nsite in Kondo model !``

*  ``ERROR ! (nelec_cond + nelec_loc + 2 * Sz) % 2 != 0 in Kondo model !``

*  ``ERROR ! nelec_cond <= nsite / 2 && 2 * |Sz| > nelec_cond + nelec_loc ...``

*  ``ERROR ! nelec_cond > nsite / 2 && abs(Sz2) > nsite / 2 * 3 - nelec...``

   In the calculation of the canonical ensemble, there are some
   irrelevant combinations of the number of electrons, the number of
   sites, and the total spin moment ( the number of electrons is larger
   twice than the number of sites); If these situations are detected,
   the program will stop.

*  | ``Check !``\ *keyword* ``is SPECIFIED but will NOT be USED.``
   | ``Please COMMENT-OUT this line``
   | ``or check this input is REALLY APPROPRIATE for your purpose !``

   Because an unnecessary parameter is specified, the program suggests
   checking the input file. If that parameter is actually unnecessary,
   please delete or comment out this line.

*  ``ERROR !``\ *keyword*\ ``is NOT specified !``

   The program stops because a prerequisite keyword is not specified.

*  *keyword* ``=`` *value* ``######  DEFAULT VALUE IS USED  ######``

   This is not an error message. The program states that the default
   value is used because this keyword is not specified.
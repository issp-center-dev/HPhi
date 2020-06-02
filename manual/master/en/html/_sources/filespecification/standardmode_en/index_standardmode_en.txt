.. _Ch:HowToStandard:

.. highlight:: none

Input files for *Standard* mode
===============================

An example of an input file for the standard mode is as follows:
::

    W = 2
    L = 4
    model = "spin"
    method = "Lanczos"

    lattice = "triangular lattice"
    //mu = 1.0
    // t = -1.0
    // t' = -0.5
    // U = 8.0
    //V = 4.0
    //V'=2.0
    J = -1.0
    J'=-0.5
    // nelec = 8
    2Sz = 0

| **Basic rules for input files**

*  In each line, there is a set of a keyword (before an “\ ``=``") and a
   parameter(after an “\ ``=``"); they are separated by “\ ``=``".

*  You can describe keywords in a random order.

*  Empty lines and lines beginning with a “``//``”(comment outs) are
   skipped.

*  Upper- and lowercase are not distinguished. Double quotes and blanks
   are ignored.

*  | There are three types of parameters.
   | 1. Parameters that must be specified (if not, :math:`{\mathcal H}\Phi` will
     stop with error messages),
   | 2. Parameters that it is not necessary to specified (if not
     specified, default values are used),
   | 3. Parameters that must not be specified (if specified,
     :math:`{\mathcal H}\Phi` will stop with error messages).
   | An example of type 3 is the transfer :math:`t` parameter for the
     Heisenberg spin system. If you choose “model=spin", you should not
     specify “\ :math:`t`".

We explain each keyword as follows.

.. toctree::
   :maxdepth: 1

   Parameters_for_the_type_of_calculation_en
   Parameters_for_the_lattice_en
   Parameters_for_conserved_quantities_en
   Parameters_for_the_Hamiltonian_en
   Parameters_for_the_numerical_condition
   Parameters_for_the_dynamical_Greens_function        
   Parameters_for_time-evolution




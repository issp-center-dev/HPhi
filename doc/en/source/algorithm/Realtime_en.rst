.. highlight:: none

Real time evolution method
==========================

In :math:`{\mathcal H}\Phi`, real time evolution calculation is done by using the following relation

.. math::

 |\Phi (t_n)\rangle = \exp (-i {\cal H}  \Delta t_n)|\Phi (t_{n-1})\rangle,

where :math:`|\Phi(t_0)\rangle` is an initial wave function and :math:`t_n = \sum_{j=1}^n  \Delta t_j`.
In calculation, we approximate :math:`\exp (-i {\cal H}  \Delta t_n)` as

.. math::

 \exp (-i {\cal H}  \Delta t_n) =\sum_{l=0}^m \frac{1}{l!}(-i {\cal H}  \Delta t_n)^l .

Here, the cut-off integer :math:`m` is set by ``ExpandCoef`` in ``ModPara``.
Each step is normalized; ``Norm.dat`` reports the norm before that
normalization. Check this value and convergence with increasing
``ExpandCoef`` to assess polynomial truncation. Energy conservation is an
additional check only for a time-independent Hamiltonian.

Time-dependent Hamiltonians
---------------------------

Each step uses one selected Hamiltonian for **all** Taylor powers,
including the first term. The time conventions are:

* ``TEOneBody`` / ``TETwoBody``: row 0 records the input state without
  propagation. Row :math:`n>0` propagates over
  :math:`\Delta t_n=t_n-t_{n-1}` with :math:`H(t_n)` (right endpoint),
  then reports the state and observables at :math:`t_n`.
  Equal consecutive times do not propagate the state.
* Raw-basis expert ``Laser``: step :math:`n` uses
  :math:`H(T_{\mathrm{init}}+n\Delta t)` and propagates for ``TimeSlice``,
  including step 0. Its output time is
  :math:`T_{\mathrm{init}}+(n+1)\Delta t`; energy uses the Hamiltonian
  selected for that step.
* Symmetry-sector ``Laser``: row 0 records the input state at ``Tinit``.
  Subsequent rows use the right-endpoint convention, as described in the
  ``TransSym`` file specification.

Standard-mode ``PumpType`` inputs that generate ``TEOneBody`` or
``TETwoBody`` follow the time-table convention, including laser pumps.
The raw expert ``Laser`` convention is distinct.

Increasing ``ExpandCoef`` improves the approximation of a frozen-H
exponential. It does not raise the order of the time discretization:
endpoint sampling is generally first order for smooth time-dependent H.
Compare successively smaller time steps at a fixed final time as well.
For example, :math:`H(t)=(1+t)K` has the exact propagator
:math:`\exp[-iK(T+T^2/2)]`; the right-endpoint error is proportional to
the step size even with a sufficiently high Taylor order.

Compatibility of earlier TE results
--------------------------------------

Releases 3.0.0 through 3.7.2 reused, in raw TE, the previous step's
:math:`H\psi` for the first Taylor term while applying the current H to
higher powers. A step that started without a prepared :math:`H\psi` did
not propagate the state. This happened in the first step of a new raw expert
``Laser`` run, and in the first step after resuming a raw TE run with
``ReStart`` = 2 or 3, even for a static Hamiltonian.
The correction changes time-dependent trajectories and resumed static
trajectories. Uninterrupted static time-table propagation is unchanged.
The correction preserves the time conventions described above.

The checkpoint format is unchanged. A resumed calculation now builds
:math:`H\psi` from the selected Hamiltonian and the saved state.
Reading a checkpoint produced by an earlier binary cannot undo errors
already accumulated in that state. Recompute affected trajectories from
their original initial state when corrected results are required.

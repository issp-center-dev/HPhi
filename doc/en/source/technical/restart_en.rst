Restart and state import
========================

HPhi has three independent mechanisms. ``ReStart`` continues a solver from
its saved algorithmic state. ``InputEigenVec=1`` imports a state for CG
evaluation or as the first state of a new TimeEvolution clock. ``CalcSpec``
4/5 resumes the Lanczos spectrum sequence through ``TMComponents`` and
``recalcvec``. Files from one mechanism must not be substituted for another.

Solver continuation
-------------------

``ReStart`` has the following values in ``CalcMod``: 0 disables continuation,
1 starts a new run and writes checkpoints, 2 reads and updates checkpoints,
and 3 reads checkpoints without rewriting them. The v3 manifest and vector
format is used by CG, mTPQ, cTPQ and TimeEvolution in raw and ``TransSym``
bases. This release requires the saved and current MPI process counts to be
equal. Sector continuation also requires the same basis layout. Rank
redistribution and layout changes are rejected.

``Lanczos_max`` is an absolute final bound for v3 continuation. CG resumes at
the saved iteration label; TPQ/cTPQ and TimeEvolution resume at the saved next
step. If the saved step is at or beyond a lowered bound, HPhi evaluates the
restored state, does not add time-series rows, and does not rewrite a
checkpoint. A terminal CG checkpoint is rechecked with the current
``LanczosEps``. A stricter tolerance can resume iterations; an upper-limit
no-op must not be described as currently converged.

``CheckpointInterval`` controls periodic CG, mTPQ and cTPQ checkpoints. Zero
means end-only. TimeEvolution uses ``OutputInterval`` for continuation
checkpoints, and its time table must contain the saved next step. See
:doc:`the restart file specification <../filespecification/outputfiles_en/restart_en>`.

State import
------------

State import validates the model, basis identity, sector and phase, payload
length, norm and digest. Source method, state, step, time and solver scalars
are provenance; they do not resume the source solver. CG evaluates an imported
state. TimeEvolution starts at step zero on its current time grid. v3 imports
currently require the same process count, and sector imports require the same
layout. The output-file reference gives the required copy procedure.

Legacy and spectrum behavior
----------------------------

Legacy raw ``tmpvec_set*_rank_*.dat`` remains readable by CG and TPQ/cTPQ.
Legacy CG keeps its historical additional-iteration meaning: iteration labels
restart at one and the old ``Lanczos_Step`` prefix is not recovered. Legacy
TPQ/cTPQ and raw TimeEvolution recover output rows below their validated next
step when possible. Old sector v1/v2 files are state-import formats, not solver
continuation checkpoints.

Ground-state Lanczos and ``CalcSpec`` 4/5 retain their existing
``TMComponents``/``recalcvec`` behavior and iteration meaning. They do not use
the v3 solver manifest.

Output recovery
---------------

Before appending, HPhi validates every required checkpoint, removes complete
rows at or after the recovery boundary, drops an incomplete final line, and
rejects malformed interior rows. A missing history file is created and only
new rows can then be written; a missing prefix cannot be reconstructed.
Per-set files and CG/TimeEvolution histories reproduce uninterrupted row
order. Shared TPQ aggregate files follow execution order, so compare them by
sorted ``(set, step, observable)`` keys when several sets resume.

The ``CalcTimer.dat`` rows ``Restart write`` (7100), ``Restart read`` (7110)
and ``Restart recover outputs`` (7120) measure v3 publication, continuation
vector loading, and time-series validation/rewrite respectively. They are
instrumentation totals, not a controlled storage benchmark.

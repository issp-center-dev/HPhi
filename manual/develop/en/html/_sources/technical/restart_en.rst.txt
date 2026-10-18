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
bases. The saved and current MPI process counts may differ. HPhi reads the
saved contiguous global intervals and redistributes their intersections to
the current ownership. Sector v3 continuation also permits metadata layout
changes between ``replicated`` (0) and ``distributed`` (1). The global vector
order remains unchanged even though rank ownership changes.

``Lanczos_max`` is an absolute final bound for v3 continuation. CG resumes at
the saved iteration label; TPQ/cTPQ and TimeEvolution resume at the saved next
step. If the bound is below the saved CG iteration label, or at or below
the saved TPQ/cTPQ/TimeEvolution next step, HPhi evaluates the restored state,
does not add time-series rows, and does not rewrite a checkpoint. CG includes
the iteration labeled ``Lanczos_max``, whereas TPQ/cTPQ/TimeEvolution exclude
that step. A terminal CG checkpoint is rechecked with the current
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
are redistributed automatically when the saved and current process counts
differ. Sector v3 and old sector v1/v2 imports also permit metadata layout
0/1 changes. For v1/v2, HPhi verifies the payload checksum using the saved
partition before redistribution and uses the partition-independent sector
digest to validate the basis set. The basis builder's canonical
representative-state sort fixes its global order, so an equal set has the
same order. The output-file reference gives the required copy procedure.

Changing the MPI process count
------------------------------

Keep every source-rank file, including a file for an empty owner. For v3
continuation, copy the manifest and every vector file named by it, then run
with the desired process count; no extra parameter is needed. v3 continuation
uses the count in each manifest row and ignores ``RestartSourceRanks``. This
also permits an interrupted TPQ/cTPQ migration to resume when some samples
are v3 files written by the new process count and other samples are legacy
files: the parameter applies only to those remaining legacy samples. It can
be removed after all samples have been republished as v3.

For v3 or sector v1/v2 state import, copy every source file under the import
prefix and run with the desired process count. Their headers provide the
source count. ``RestartSourceRanks`` is optional for state import and acts
only as an additional equality check.

For legacy raw import or continuation, copy every saved per-rank file and add
``RestartSourceRanks=<saved process count>`` to ``ModPara`` before changing
the process count. HPhi never guesses this value from file names, because a
missing empty trailing rank would otherwise be indistinguishable from a
smaller complete file set.

Legacy and spectrum behavior
----------------------------

Legacy raw ``tmpvec_set*_rank_*.dat`` remains readable by CG and TPQ/cTPQ.
Because this format stores neither its source rank count nor identity words,
rank-changing reads require ``RestartSourceRanks``. Keep every saved rank
file, including an empty trailing rank; HPhi checks all native headers, file
lengths, steps and the global dimension before redistribution, but cannot
validate model or sector identity from the files.
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

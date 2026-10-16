"""CheckpointInterval / RestartSourceRanks parsing and combination rules."""
import restart_common as rc
ROOT = rc.fresh_dir("restart_input_validation")
CG = 'L=8\nmodel="Spin"\nmethod="CG"\nlattice="chain"\nJ=1\n2Sz=0\nexct=1\nLanczos_max=20\noutputmode="None"\n'
case = ROOT / "cg"; rc.generate(case, CG)
rc.set_key(case, "modpara.def", "CheckpointInterval", 5)
rc.run(case, "interval_without_restart", fail="CheckpointInterval requires ReStart=1 or 2", mpi=False)
rc.set_key(case, "calcmod.def", "ReStart", 3)
rc.run(case, "interval_with_in_only", fail="CheckpointInterval requires ReStart=1 or 2", mpi=False)
rc.set_key(case, "modpara.def", "CheckpointInterval", -1)
rc.set_key(case, "calcmod.def", "ReStart", 1)
rc.run(case, "negative_interval", fail="CheckpointInterval must be a non-negative integer", mpi=False)
rc.set_key(case, "modpara.def", "CheckpointInterval", 0)
rc.set_key(case, "modpara.def", "RestartSourceRanks", -2)
rc.run(case, "negative_ranks", fail="RestartSourceRanks must be a non-negative integer", mpi=False)

for key, value, label in (
    ("CheckpointInterval", 1.5, "fractional_interval"),
    ("CheckpointInterval", "nan", "nonfinite_interval"),
    ("CheckpointInterval", 2147483647, "int_max_interval"),
    ("RestartSourceRanks", 1.5, "fractional_ranks"),
    ("RestartSourceRanks", "nan", "nonfinite_ranks"),
    ("RestartSourceRanks", 2147483647, "int_max_ranks"),
):
    rc.set_key(case, "modpara.def", "CheckpointInterval", 0)
    rc.set_key(case, "modpara.def", "RestartSourceRanks", 0)
    rc.set_key(case, "modpara.def", key, value)
    rc.run(case, label, fail=key + " must be a non-negative integer", mpi=False)

# Explicit zero is valid for both keys and preserves the disabled defaults.
rc.set_key(case, "calcmod.def", "ReStart", 0)
rc.set_key(case, "modpara.def", "CheckpointInterval", 0)
rc.set_key(case, "modpara.def", "RestartSourceRanks", 0)
rc.set_key(case, "modpara.def", "Lanczos_max", 100)
rc.run(case, "explicit_zeros", mpi=False)

# Standard mode writes the keys only when given.
plain = ROOT / "plain"; rc.generate(plain, CG)
assert "CheckpointInterval" not in (plain / "modpara.def").read_text()
assert "RestartSourceRanks" not in (plain / "modpara.def").read_text()
given = ROOT / "given"; rc.generate(given, CG + "CheckpointInterval=3\nRestartSourceRanks=4\n")
text = (given / "modpara.def").read_text()
assert "CheckpointInterval 3" in " ".join(text.split()) and "RestartSourceRanks 4" in " ".join(text.split()), text
# TE rejects a positive CheckpointInterval and ignores zero (TE run itself is exercised in Task 14).
te = ROOT / "te"; rc.generate(te, CG)
rc.set_key(te, "calcmod.def", "CalcType", 4)
rc.set_key(te, "calcmod.def", "ReStart", 1)
rc.set_key(te, "modpara.def", "CheckpointInterval", 2)
rc.run(te, "te_interval", fail="CheckpointInterval is not used by TimeEvolution", mpi=False)
print("ok")

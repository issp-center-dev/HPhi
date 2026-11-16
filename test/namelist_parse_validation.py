"""Exercise namelist parsing through Expert mode, before any calculation."""
from pathlib import Path
import os
import subprocess
import sys

HPHI = str(Path(sys.argv[1]).resolve())
ROOT = Path("namelist_parse_validation")
ROOT.mkdir(exist_ok=True)
# This intentional downstream error proves that GetFileName completed and
# selected the requested CalcMod file. No model setup or MPI launch is needed.
SENTINEL = "CalcType: 0: Lanczos Method"
BASE = "CalcMod calcmod.def\nModPara unused.def\nLocSpin unused.def\n"
PAIR = "keyword and filename must be set as a pair"
LONG_FILE = "namelist filename exceeds 199 bytes"
LONG_LINE = "namelist line exceeds 255 bytes"


def check(name, namelist, expected, calcmod="calcmod.def"):
    path = ROOT / name
    path.mkdir(exist_ok=True)
    (path / "namelist.def").write_bytes(namelist.encode("ascii"))
    (path / calcmod).write_text("CalcType 99\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    # Intentionally ignore MPIRUN: these are serial parser/entry-guard tests,
    # including when CTest is invoked from the np16 baseline job.
    with (path / "run.log").open("w") as log:
        result = subprocess.run(
            [HPHI, "-e", "namelist.def"], cwd=path, env=env,
            stdout=log, stderr=subprocess.STDOUT, timeout=30)
    output = (path / "run.log").read_text(errors="replace")
    if result.returncode <= 0 or expected not in output:
        raise AssertionError(
            f"{name}: rc={result.returncode}, expected {expected!r}\n{output}")
    if expected != SENTINEL and SENTINEL in output:
        raise AssertionError(f"{name}: reached CalcMod after malformed input\n{output}")
    print(f"PASS {name}")


check("plain", BASE, SENTINEL)
check("reordered_case", "locspin unused.def\nMODPARA unused.def\ncalcmod calcmod.def\n",
      SENTINEL)
check("blank_comments", "# header\n\n \t\r\n  # indented\n" +
      BASE.replace("\n", "\n\t\n # between pairs\n") + "# eof", SENTINEL)
check("whitespace", " \tCalcMod\tcalcmod.def\r\nModPara\vunused.def\n"
      "LocSpin\funused.def", SENTINEL)
check("no_final_newline", BASE.rstrip("\n"), SENTINEL)
check("trailing_words", BASE.replace("calcmod.def", "calcmod.def # inline comment")
      .replace("LocSpin unused.def", "LocSpin unused.def ignored"), SENTINEL)
check("missing_first_filename", "CalcMod\n" + BASE, PAIR)
check("missing_later_filename", BASE + "Trans\n", PAIR)
check("missing_filename_spaces", BASE + "Trans \t\r\n", PAIR)
check("missing_filename_eof", BASE + "Trans", PAIR)
check("filename_without_keyword", "calcmod.def\n" + BASE, PAIR)
check("unknown_keyword", BASE + "Unknown unused.def\n", "Wrong keywords 'Unknown'")
check("long_keyword", "K" * 200 + " unused.def\n", "Wrong keywords")
check("duplicate_keyword", BASE + "CalcMod other.def\n", "Same keywords exist")
check("duplicate_case", BASE + "calcmod other.def\n", "Same keywords exist")
check("empty", "", "Need to make a def file for CalcMod")
check("comments_only", "# a\n \t\n # b", "Need to make a def file for CalcMod")
long_name = "c" * 199
check("max_filename", BASE.replace("calcmod.def", long_name), SENTINEL, long_name)
check("oversize_filename", BASE.replace("calcmod.def", "c" * 200), LONG_FILE)
check("oversize_last_filename", BASE + "Trans " + "t" * 200, LONG_FILE)
# A valid-looking prefix must not hide an overlong physical line.
line = "CalcMod calcmod.def"
for suffix, newline in (("eof", ""), ("lf", "\n"), ("crlf", "\r\n")):
    prefix = "ModPara unused.def\nLocSpin unused.def\n"
    check("max_line_" + suffix, prefix + line.ljust(255) + newline, SENTINEL)
    check("oversize_line_" + suffix, prefix + line.ljust(256) + newline, LONG_LINE)
check("oversize_comment", "#" + "c" * 255 + "\n" + BASE, LONG_LINE)
print("Namelist parser: 27 cases passed.")

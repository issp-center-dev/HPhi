#!/usr/bin/env python3
# HPhi ExpecLocal call-inventory tooling.
#
# Strip C comments from a source file using a small state machine (not a
# regex) so that "//" or "/*" occurring inside string or character literals
# is never mistaken for the start of a comment. Used by
# test/check_expec_local_calls.sh (and by the Step 1 inventory generation in
# docs/superpowers/specs/2026-07-11-expec-call-inventory.md) so that grep-based
# scanning for MPI call tokens does not see commented-out code.
#
# gcc -fpreprocessed is not available on this project's macOS/Apple clang
# toolchain, hence this standalone stripper instead of relying on the
# preprocessor to do comment removal.
import sys

src = open(sys.argv[1]).read()
out, i, n, st = [], 0, len(src), 'code'   # st: code|blk|line|str|chr
while i < n:
    c = src[i]; c2 = src[i:i+2]
    if st == 'code':
        if c2 == '/*': st = 'blk'; i += 2; continue
        if c2 == '//': st = 'line'; i += 2; continue
        if c == '"': st = 'str'
        elif c == "'": st = 'chr'
        out.append(c); i += 1
    elif st == 'blk':
        if c == '\n': out.append(c)
        if c2 == '*/': st = 'code'; i += 2
        else: i += 1
    elif st == 'line':
        if c == '\n': st = 'code'; out.append(c)
        i += 1
    else:  # str / chr
        if c == '\\': out.append(src[i:i+2]); i += 2; continue
        if (st == 'str' and c == '"') or (st == 'chr' and c == "'"): st = 'code'
        out.append(c); i += 1
res = ''.join(out)
assert res.strip(), f"stripper produced empty output for {sys.argv[1]}"
sys.stdout.write(res)

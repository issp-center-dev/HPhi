/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/**
 * @file green_row_format.h
 *
 * @brief Single source of truth for the one-body and two-body Green's-
 * function data-row format strings, shared by src/expec_cisajs.c /
 * src/expec_cisajscktaltdc.c (ExpecMode 1) and src/expec_trace.c (ExpecMode
 * 2 trace-kernel output phases, phase 3b Tasks 3/4).
 *
 * Every expec_cisajs_* writer prints one row per (i,sigma,j,sigma) pair as
 * `" %4lu %4lu %4lu %4lu %.10lf %.10lf\n"` with args (isite1-1, sigma1,
 * isite2-1, sigma2, Re(value), Im(value)), immediately preceded by a
 * GreenOutputWriteIndexPrefix(fp, X) call. The trace kernel's output phase
 * must produce byte-identical rows for the same values, so both call sites
 * use this one macro instead of two independently-typed literals that could
 * drift apart.
 *
 * The two-body row is the same idea (8 index fields instead of 4) but
 * src/expec_cisajscktaltdc.c's PRE-EXISTING call sites are not all
 * byte-identical to each other: every site uses "%4ld" (not "%4lu") for the
 * index fields, but they split into two groups by whether a single trailing
 * space precedes the final "\n". GREEN_TWOBODY_ROW_FORMAT (no trailing
 * space) is what expec_cisajscktalt_HubbardGC's row and
 * expec_cisajscktalt_Hubbard's normally-computed row use.
 * GREEN_TWOBODY_ROW_FORMAT_SP (one trailing space before "\n") is what
 * every other reachable two-body row in that file uses: canonical Hubbard's
 * Sz-conserved-violation 0.0 shortcut row, and every row
 * expec_cisajscktalt_SpinHalf/SpinGCHalf write (both their Rearray-
 * irregular 0.0 row and their normally-computed row). This split predates
 * phase 3b and is preserved verbatim rather than "fixed", since the trace
 * kernel's job is byte-identical output, not reformatting Mode 1.
 */
#pragma once

#define GREEN_ONEBODY_ROW_FORMAT " %4lu %4lu %4lu %4lu %.10lf %.10lf\n"

#define GREEN_TWOBODY_ROW_FORMAT " %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf\n"
#define GREEN_TWOBODY_ROW_FORMAT_SP " %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %.10lf %.10lf \n"

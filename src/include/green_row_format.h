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
 * @brief Single source of truth for the one-body Green's-function data-row
 * format string, shared by src/expec_cisajs.c (ExpecMode 1) and
 * src/expec_trace.c (ExpecMode 2 trace-kernel one-body output phase, phase
 * 3b Task 3).
 *
 * Every expec_cisajs_* writer prints one row per (i,sigma,j,sigma) pair as
 * `" %4lu %4lu %4lu %4lu %.10lf %.10lf\n"` with args (isite1-1, sigma1,
 * isite2-1, sigma2, Re(value), Im(value)), immediately preceded by a
 * GreenOutputWriteIndexPrefix(fp, X) call. The trace kernel's output phase
 * must produce byte-identical rows for the same values, so both call sites
 * use this one macro instead of two independently-typed literals that could
 * drift apart.
 */
#pragma once

#define GREEN_ONEBODY_ROW_FORMAT " %4lu %4lu %4lu %4lu %.10lf %.10lf\n"

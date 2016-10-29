/*
ISSP Math Library - A library for solving linear systems in materials science
Copyright (C) 2016 Mitsuaki Kawamura

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

For more details, See ÅeCOPYING.LESSERÅf in the root directory of this library.
*/

#include <complex.h>
#pragma once

#ifdef INTEL
#define komega_bicg_init komega_bicg_mp_komega_bicg_init_
#define komega_bicg_restart komega_bicg_mp_komega_bicg_restart_
#define komega_bicg_update komega_bicg_mp_komega_bicg_update_
#define komega_bicg_getcoef komega_bicg_mp_komega_bicg_getcoef_
#define komega_bicg_getvec komega_bicg_mp_komega_bicg_getvec_
#define komega_bicg_finalize komega_bicg_mp_komega_bicg_finalize_
#else
#define komega_bicg_init __komega_bicg_MOD_komega_bicg_init
#define komega_bicg_restart __komega_bicg_MOD_komega_bicg_restart
#define komega_bicg_update __komega_bicg_MOD_komega_bicg_update
#define komega_bicg_getcoef __komega_bicg_MOD_komega_bicg_getcoef
#define komega_bicg_getvec __komega_bicg_MOD_komega_bicg_getvec
#define komega_bicg_finalize __komega_bicg_MOD_komega_bicg_finalize
#endif

void komega_bicg_init(int *ndim, int *nl, int *nz, double complex *x, double complex *z, int *itermax, double *threshold);
void komega_bicg_restart(int *ndim, int *nl, int *nz, double complex *x, double complex *z, int *itermax, double *threshold, int *status,
  int *iter_old, double complex *v2, double complex *v12, double complex *v4, double complex *v14, double complex *alpha_save, double complex *beta_save, double complex *z_seed, double complex *r_l_save);
void komega_bicg_update(double complex *v12, double complex *v2, double complex *v14, double complex *v4, double complex *x, double complex *r_l, int *status);
void komega_bicg_getcoef(double complex *alpha_save, double complex *beta_save, double complex *z_seed, double complex *r_l_save);
void komega_bicg_getvec(double complex *r_old, double complex *r_tilde_old);
void komega_bicg_finalize();

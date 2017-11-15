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

For more details, See ‘COPYING.LESSER’ in the root directory of this library.
*/

#include <complex.h>
#pragma once

void komega_bicg_init(int *ndim, int *nl, int *nz, double complex *x,
                      double complex *z, int *itermax, double *threshold, int *comm);
void komega_cocg_init(int *ndim, int *nl, int *nz, double complex *x,
                      double complex *z, int *itermax, double *threshold, int *comm);
void komega_cg_c_init(int *ndim, int *nl, int *nz, double complex *x,
                      double *z, int *itermax, double *threshold, int *comm);
void komega_cg_r_init(int *ndim, int *nl, int *nz, double *x,
                      double *z, int *itermax, double *threshold, int *comm);

void komega_bicg_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double complex *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, double complex *v4, double complex *v14,
                         double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save, int *comm);
void komega_cocg_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double complex *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, 
                         double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save, int *comm);
void komega_cg_c_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, 
                         double *alpha_save, double *beta_save,
                         double *z_seed, double complex *r_l_save, int *comm);
void komega_cg_r_restart(int *ndim, int *nl, int *nz, double *x,
                         double *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double *v2,
                         double *v12, 
                         double *alpha_save, double *beta_save,
                         double *z_seed, double *r_l_save, int *comm);

void komega_bicg_update(double complex *v12, double complex *v2,
                        double complex *v14, double complex *v4,
                        double complex *x, double complex *r_l, int *status);
void komega_cocg_update(double complex *v12, double complex *v2,
                        double complex *x, double complex *r_l, int *status);
void komega_cg_c_update(double complex *v12, double complex *v2,
                        double complex *x, double complex *r_l, int *status);
void komega_cg_r_update(double *v12, double *v2,
                        double *x, double *r_l, int *status);


void komega_bicg_getcoef(double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save);
void komega_cocg_getcoef(double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save);
void komega_cg_c_getcoef(double *alpha_save, double *beta_save,
                         double *z_seed, double complex *r_l_save);
void komega_cg_r_getcoef(double *alpha_save, double *beta_save,
                         double *z_seed, double *r_l_save);

void komega_bicg_getvec(double complex *r_old, double complex *r_tilde_old);
void komega_cocg_getvec(double complex *r_old);
void komega_cg_c_getvec(double complex *r_old);
void komega_cg_r_getvec(double *r_old);

void komega_bicg_getresidual(double *res);
void komega_cocg_getresidual(double *res);
void komega_cg_c_getresidual(double *res);
void komega_cg_r_getresidual(double *res);

void komega_bicg_finalize();
void komega_cocg_finalize();
void komega_cg_r_finalize();
void komega_cg_c_finalize();

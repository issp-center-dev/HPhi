/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <complex.h>
#include <stdio.h>

void StdFace_intr(struct StdIntList *StdI, double complex intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

void StdFace_Hopping(struct StdIntList *StdI, double complex trans0,
  int isite, int jsite);
void StdFace_MagField(struct StdIntList *StdI, int S2, double h, double Gamma, int isite);

void StdFace_Coulomb(struct StdIntList *StdI, double V, int isite, int jsite);
void StdFace_GeneralJ(struct StdIntList *StdI, double J[3][3],
  int Si2, int Sj2, int isite, int jsite);

void StdFace_PrintVal_d(char* valname, double *val, double val0);
void StdFace_PrintVal_dd(char* valname, double *val, double val0, double val1);
void StdFace_PrintVal_c(char* valname, double complex *val, double complex val0);
void StdFace_PrintVal_i(char* valname, int *val, int val0);

void StdFace_NotUsed_d(char* valname, double val);
void StdFace_NotUsed_i(char* valname, int val);
void StdFace_NotUsed_c(char* valname, double complex val);
void StdFace_NotUsed_J(char* valname, double JAll, double J[3][3]);

void StdFace_RequiredVal_i(char* valname, int val);
void StdFace_InputSpinNN(struct StdIntList *StdI, double J0[3][3],
  double J0All, char *J0name);
void StdFace_InputSpin(struct StdIntList *StdI, double Jp[3][3],
  double JpAll, char *Jpname);
void StdFace_InputCoulombV(struct StdIntList *StdI, double *V0, char *V0name);
void StdFace_InputHopp(struct StdIntList *StdI, double complex *t0, char *t0name);

void StdFace_InitSite2D(struct StdIntList *StdI, FILE *fp);
void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp,
  int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC,
  int *isite, int *jsite, int connect);
void StdFace_InitSite1D(struct StdIntList *StdI);
void StdFace_PrintGeometry(struct StdIntList *StdI);

void StdFace_Tetragonal(struct StdIntList *StdI, char *model);
void StdFace_Chain(struct StdIntList *StdI, char *model);
void StdFace_Ladder(struct StdIntList *StdI, char *model);
void StdFace_Triangular(struct StdIntList *StdI, char *model);
void StdFace_Honeycomb(struct StdIntList *StdI, char *model);
void StdFace_Kagome(struct StdIntList *StdI, char *model);

void StdFace_Chain_Boost(struct StdIntList *StdI);
void StdFace_Ladder_Boost(struct StdIntList *StdI);
void StdFace_Honeycomb_Boost(struct StdIntList *StdI);
void StdFace_Kagome_Boost(struct StdIntList *StdI);

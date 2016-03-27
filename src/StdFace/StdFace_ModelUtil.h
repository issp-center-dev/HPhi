/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima

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

void StdFace_intr(struct StdIntList *StdI, _Dcomplex intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

void StdFace_Hopping(struct StdIntList *StdI, _Dcomplex trans0,
  int isite, int jsite, char *mode);
void StdFace_MagField(struct StdIntList *StdI, int S2, double h, double Gamma, int isite);

void StdFace_Coulomb(struct StdIntList *StdI, double V, int isite, int jsite);
void StdFace_GeneralJ(struct StdIntList *StdI, double *J,
  int Si2, int Sj2, int isite, int jsite);

void StdFace_PrintVal_d(char* valname, double *val, double val0);
void StdFace_PrintVal_c(char* valname, _Dcomplex *val, _Dcomplex val0);
void StdFace_PrintVal_i(char* valname, int *val, int val0);
void StdFace_NotUsed_d(char* valname, double val);
void StdFace_NotUsed_i(char* valname, int val);
void StdFace_NotUsed_c(char* valname, _Dcomplex val);
void StdFace_RequiredVal_i(char* valname, int val);

void StdFace_InitSite2D(struct StdIntList *StdI, FILE *fp);
void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp,
  int iW, int iL, int jW, int jL,
  double xiW, double xiL, double xjW, double xjL,
  int isite, int *jsite, int connect, int NsiteUC, int isiteUC, char *model);

void StdFace_SquareLattice(struct StdIntList *StdI, char *model);
void StdFace_ChainLattice(struct StdIntList *StdI, char *model);
void StdFace_ChainLattice_Boost(struct StdIntList *StdI, int Sz2, int lGC);
void StdFace_Ladder(struct StdIntList *StdI, char *model);
void StdFace_Ladder_Boost(struct StdIntList *StdI, int Sz2, int lGC);
void StdFace_TriangularLattice(struct StdIntList *StdI, char *model);

void FermionHubbard_HoneycombLattice(struct StdIntList *StdI, int nelec, int lGC);
void KondoLattice_HoneycombLattice(struct StdIntList *StdI, int nelec, int lGC);
void Spin_HoneycombLattice(struct StdIntList *StdI, int Sz2, int lGC);
void Spin_HoneycombLattice_Boost(struct StdIntList *StdI, int Sz2, int lGC);

void Spin_Kagome_Boost(struct StdIntList *StdI, int Sz2, int lGC);

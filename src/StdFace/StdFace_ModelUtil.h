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
void StdFace_trans(int *ktrans, double trans0,
  int isite, int ispin, int jsite, int jspin);

void StdFace_intr(int *kintr, double intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

void StdFace_exchange(int *kintr, double Jexc, int isite, int jsite);
void StdFace_Kitaev(int *kintr, double Jflip, int isite, int jsite);
void StdFace_SzSz(int *kintr, double Jexc, int isite, int jsite);
void StdFace_Coulomb(int *kintr, double V, int isite, int jsite);

void StdFace_PrintVal_d(char* valname, double *val, double val0);
void StdFace_PrintVal_i(char* valname, int *val, int val0);
void StdFace_NotUsed_d(char* valname, double val);
void StdFace_NotUsed_i(char* valname, int val);
void StdFace_RequiredVal_i(char* valname, int val);

void FermionHubbard_SquareLattice(int nelec, int lGC);
void KondoLattice_SquareLattice(int nelec, int lGC);
void Spin_SquareLattice();

void FermionHubbard_ChainLattice(int nelec, int lGC);
void KondoLattice_ChainLattice(int nelec, int lGC);
void Spin_ChainLattice();

void FermionHubbard_TriangularLattice(int nelec, int lGC);
void KondoLattice_TriangularLattice(int nelec, int lGC);
void Spin_TriangularLattice();

void FermionHubbard_HoneycombLattice(int nelec, int lGC);
void KondoLattice_HoneycombLattice(int nelec, int lGC);
void Spin_HoneycombLattice();

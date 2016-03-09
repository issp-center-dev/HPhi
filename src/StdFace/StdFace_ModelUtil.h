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
void StdFace_trans(struct StdIntList *StdI, double trans0,
  int isite, int ispin, int jsite, int jspin);

void StdFace_intr(struct StdIntList *StdI, double intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

void StdFace_MagLong(struct StdIntList *StdI, int S2, double Mag, int isite);
void StdFace_MagTrans(struct StdIntList *StdI, int S2, double Mag, int isite);

void StdFace_exchange(struct StdIntList *StdI, int Si2, int Sj2, double Jexc, int isite, int jsite);
void StdFace_Kitaev(struct StdIntList *StdI, int Si2, int Sj2, double Jflip, int isite, int jsite);
void StdFace_SzSz(struct StdIntList *StdIr, int Si2, int Sj2, double Jexc, int isite, int jsite);
void StdFace_Coulomb(struct StdIntList *StdI, double V, int isite, int jsite);

void StdFace_PrintVal_d(char* valname, double *val, double val0);
void StdFace_PrintVal_i(char* valname, int *val, int val0);
void StdFace_NotUsed_d(char* valname, double val);
void StdFace_NotUsed_i(char* valname, int val);
void StdFace_RequiredVal_i(char* valname, int val);

void FermionHubbard_SquareLattice(struct StdIntList *StdI, int nelec, int lGC);
void KondoLattice_SquareLattice(struct StdIntList *StdI, int nelec, int lGC);
void Spin_SquareLattice(struct StdIntList *StdI, int Sz2, int lGC);

void FermionHubbard_ChainLattice(struct StdIntList *StdI, int nelec, int lGC);
void KondoLattice_ChainLattice(struct StdIntList *StdI, int nelec, int lGC);
void Spin_ChainLattice(struct StdIntList *StdI, int Sz2, int lGC);
void Spin_ChainLattice_Boost(struct StdIntList *StdI, int Sz2, int lGC);

void FermionHubbard_TriangularLattice(struct StdIntList *StdI, int nelec, int lGC);
void KondoLattice_TriangularLattice(struct StdIntList *StdI, int nelec, int lGC);
void Spin_TriangularLattice(struct StdIntList *StdI, int Sz2, int lGC);

void FermionHubbard_HoneycombLattice(struct StdIntList *StdI, int nelec, int lGC);
void KondoLattice_HoneycombLattice(struct StdIntList *StdI, int nelec, int lGC);
void Spin_HoneycombLattice(struct StdIntList *StdI, int Sz2, int lGC);
void Spin_HoneycombLattice_Boost(struct StdIntList *StdI, int Sz2, int lGC);

void FermionHubbard_Ladder(struct StdIntList *StdI, int nelec, int lGC);
void KondoLattice_Ladder(struct StdIntList *StdI, int nelec, int lGC);
void Spin_Ladder(struct StdIntList *StdI, int Sz2, int lGC);
void Spin_Ladder_Boost(struct StdIntList *StdI, int Sz2, int lGC);

void Spin_Kagome_Boost(struct StdIntList *StdI, int Sz2, int lGC);

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

/*
Parameters for LATTICE
*/
double a; /**< The lattice constant */
int L;
int W;
/*
Parameters for MODEL
*/
double mu;
double t;
double tp;
double tpp;
double t0;
double t1;
double t1p;
double t2;
double U;
double V;
double Vp;
double Vpp;
double V0;
double V1;
double V1p;
double V2;
double J;
double Jp;
double Jpp;
double J0;
double J1;
double J1p;
double J2;
double J2p;
/**/
double Jx;
double Jy;
double Jz;
double Jx0;
double Jy0;
double Jz0;
double Jx1;
double Jy1;
double Jz1;
double Jx2;
double Jy2;
double Jz2;
double Jxp;
double Jyp;
double Jzp;
double Jxy;
double Jxy0;
double Jxy1;
double Jxy2;
double Jxyp;
double h;
double Gamma;
double D;
double K;

int nsite;
int *locspinflag;
int ntrans;
int **transindx;
double *trans;
int nintr;
int **intrindx;
double *intr;

int LargeValue;
int S2;

int ***list_6spin_pair;
int **list_6spin_star;
int num_pivot;
int ishift_nspin;
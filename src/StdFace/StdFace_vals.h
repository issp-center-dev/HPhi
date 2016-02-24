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
static double a; /**< The lattice constant */
static int L;
static int W;
/*
Parameters for MODEL
*/
static double mu;
static double t;
static double tp;
static double tpp;
static double t0;
static double t1;
static double t1p;
static double t2;
static double U;
static double V;
static double Vp;
static double Vpp;
static double V0;
static double V1;
static double V1p;
static double V2;
static double J;
static double Jp;
static double Jpp;
static double J0;
static double J1;
static double J1p;
static double J2;
static double J2p;
/**/
static double Jx;
static double Jy;
static double Jz;
static double Jx0;
static double Jy0;
static double Jz0;
static double Jx1;
static double Jy1;
static double Jz1;
static double Jx2;
static double Jy2;
static double Jz2;
static double Jxp;
static double Jyp;
static double Jzp;
static double Jxy;
static double Jxy0;
static double Jxy1;
static double Jxy2;
static double Jxyp;
static double h;
static double Gamma;
static double D;
static double K;

static int nsite;
static int *locspinflag;
static int ntrans;
static int **transindx;
static double *trans;
static int nintr;
static int **intrindx;
static double *intr;

static int LargeValue;
static int S2;

static int ***list_6spin_pair;
static int **list_6spin_star;
static int num_pivot;
static int ishift_nspin;
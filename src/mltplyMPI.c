
#include "wrapperMPI.h"
#include "mpi.h"
#include <complex.h>
#include <stdlib.h>
#include "struct.h"
#include "global.h"
#include "mltply.h"

void child_general_hopp_MPIdouble(int itrans, struct BindStruct *X, 
  double complex *tmp_v0, double complex *tmp_v1)
{
  int bit1, bit2, mybit1, mybit2, ierr, dest, bitdiff, sgn;
  unsigned long int idim_max_buf, j;
  MPI_Status status;
  double complex trans, dam_pr0, dam_pr;

  bit1 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][0] - 2 
                            + X->Def.EDGeneralTransfer[itrans][1]];
  bit2 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][2] - 2
                            + X->Def.EDGeneralTransfer[itrans][3]];
  bitdiff = abs(bit1 - bit2);

  mybit1 = myrank & bit1;
  mybit2 = myrank & bit2;

  dest = myrank ^ (bit1 + bit2);
  SgnBit((unsigned long int)(dest & bitdiff), &sgn); // Fermion sign

  if(mybit1 == 0 && mybit2 == bit2){
    trans = (double)sgn * X->Def.EDParaGeneralTransfer[itrans];
  }
  else if(mybit1 == bit1 && mybit2 == 0) {
    trans = (double)sgn * conj(X->Def.EDParaGeneralTransfer[itrans]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, dest, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, dest, 0, MPI_COMM_WORLD, &status);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max, MPI_DOUBLE_COMPLEX, dest, 0,
    v1buf, idim_max_buf, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD, &status);

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY) { // for multply
    for (j = 0; j < idim_max_buf; j++) {
      dam_pr0 = trans * v1buf[j];
      tmp_v0[j] += dam_pr0;
      dam_pr += tmp_v1[j] * dam_pr0;
    }
  }
  else {
    for (j = 0; j < idim_max_buf; j++) 
      dam_pr += tmp_v1[j] * trans * v1buf[j];
  }

  X->Large.prdct += dam_pr;
  
}

void child_general_hopp_MPIsingle(int itrans, struct BindStruct *X,
  double complex *tmp_v0, double complex *tmp_v1)
{
  int bit1, bit2, mybit1, mybit2, ierr, dest, bitdiff, sgn;
  unsigned long int idim_max_buf, j;
  MPI_Status status;
  double complex trans, dam_pr0, dam_pr;

  bit1 = X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][0] - 2
    + X->Def.EDGeneralTransfer[itrans][1]];
  bit2 = (int)X->Def.Tpow[2 * X->Def.EDGeneralTransfer[itrans][2] - 2
    + X->Def.EDGeneralTransfer[itrans][3]];
  bitdiff = bit2 - 1;

  mybit2 = myrank & bit2;

  dest = myrank ^ bit2;
  SgnBit((unsigned long int)(dest & bitdiff), &sgn); // Fermion sign

  if (mybit1 == 0 && mybit2 == bit2) {
    trans = (double)sgn * X->Def.EDParaGeneralTransfer[itrans];
  }
  else if (mybit1 == bit1 && mybit2 == 0) {
    trans = (double)sgn * conj(X->Def.EDParaGeneralTransfer[itrans]);
  }
  else return;

  ierr = MPI_Sendrecv(&X->Check.idim_max, 1, MPI_UNSIGNED_LONG, dest, 0,
    &idim_max_buf, 1, MPI_UNSIGNED_LONG, dest, 0, MPI_COMM_WORLD, &status);
  ierr = MPI_Sendrecv(tmp_v1, X->Check.idim_max, MPI_DOUBLE_COMPLEX, dest, 0,
    v1buf, idim_max_buf, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD, &status);

  dam_pr = 0.0;
  if (X->Large.mode == M_MLTPLY) { // for multply
    for (j = 0; j < idim_max_buf; j++) {
      dam_pr0 = trans * v1buf[j];
      tmp_v0[j] += dam_pr0;
      dam_pr += tmp_v1[j] * dam_pr0;
    }
  }
  else {
    for (j = 0; j < idim_max_buf; j++)
      dam_pr += tmp_v1[j] * trans * v1buf[j];
  }

  X->Large.prdct += dam_pr;

}


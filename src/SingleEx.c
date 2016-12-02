#include "bitcalc.h"
#include "SingleEx.h"
#include "wrapperMPI.h"
#include "mfmemory.h"
#include "mltplyMPI.h"
#ifdef MPI
#include "mpi.h"
#endif

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param is1_spin
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 */
  double complex GC_Cis(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          long unsigned int is1_spin,
          double complex tmp_V,
          long unsigned int *tmp_off
  ) {

    long unsigned int list_1_j, list_1_off;
    long unsigned int ibit_tmp_1;
    long unsigned int bit;
    int sgn, ipsgn;
    double complex dmv, dam_pr;

    list_1_j = j - 1;

    ibit_tmp_1 = (list_1_j & is1_spin);
    // is1_spin >= 1
    // is1_spin = Tpow[2*isite + ispin]

    *tmp_off = 0;

    if (ibit_tmp_1 == 0) {
    // able to create an electron at the is1_spin state
      bit = list_1_j - ( list_1_j & (2*is1_spin-1) );
      SgnBit(bit, &sgn); // Fermion sign 
      ipsgn = 1;
#ifdef MPI
      SgnBit(myrank, &ipsgn); // Fermion sign 
#endif
      list_1_off = list_1_j | is1_spin; // OR
      *tmp_off = list_1_off;
      dmv = ipsgn * sgn * tmp_v1[j];
      //if (X->Large.mode == M_MLTPLY) { // for multply
      tmp_v0[list_1_off + 1] += dmv * tmp_V;
      //}
      dam_pr = dmv * conj(tmp_v1[list_1_off + 1]);
      return dam_pr;
    } else {
      return 0;
    }
  }

/**
 *
 *
 * @param j
 * @param tmp_v0
 * @param tmp_v1
 * @param X
 * @param is1_spin
 * @param tmp_V
 *
 * @return
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 */
  double complex GC_Ajt(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          long unsigned int is1_spin,
          double complex tmp_V,
          long unsigned int *tmp_off
  ) {

    long unsigned int list_1_j, list_1_off;
    long unsigned int ibit_tmp_1;
    long unsigned int bit;
    int sgn, ipsgn;
    double complex dmv, dam_pr;

    list_1_j = j - 1;

    ibit_tmp_1 = (list_1_j & is1_spin);
    // is1_spin >= 1

    *tmp_off = 0;

    if (ibit_tmp_1 == is1_spin) {
    // able to create an electron at the is1_spin state
      bit = list_1_j - ( list_1_j & (2*is1_spin-1) );
      SgnBit(bit, &sgn); // Fermion sign 
      ipsgn = 1;
#ifdef MPI
      SgnBit(myrank, &ipsgn); // Fermion sign 
#endif
      list_1_off = list_1_j ^ is1_spin;
      *tmp_off = list_1_off;
      dmv = ipsgn * sgn * tmp_v1[j];
      //if (X->Large.mode == M_MLTPLY) { // for multply
      tmp_v0[list_1_off + 1] += dmv * tmp_V;
      //}
      dam_pr = dmv * conj(tmp_v1[list_1_off + 1]);
      return dam_pr;
    } else {
      return 0;
    }
  }


/**
  *
  * Single creation/annihilation operator
  * in the inter process region for HubbardGC.
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  * @author Kazuyoshi Yoshimi (The University of Tokyo)
  * @author Youhei Yamaji (The University of Tokyo)
  */
double complex X_GC_Cis_MPI(
				       int org_isite,
				       int org_ispin,
				       double complex tmp_trans,
  double complex *tmp_v0 /**< [out] Result v0 += H v1*/,
  double complex *tmp_v1 /**< [in] v0 += H v1*/,
  unsigned long int idim_max,
  double complex *tmp_v1buf,
  long unsigned int *Tpow
  ) {
#ifdef MPI
    int mask2, state2, ierr, origin, bit2diff, Fsgn;
    unsigned long int idim_max_buf, j;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    // org_isite >= Nsite
    mask2 = (int) Tpow[2 * org_isite + org_ispin];

    origin = myrank ^ mask2; // XOR
    state2 = origin & mask2;

    //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
    //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin. 

    bit2diff = myrank - ((2*mask2-1) & myrank);

    //SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign
    SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign

    ierr = MPI_Sendrecv(&idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        tmp_v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    if (state2 == mask2) {
        trans = 0;
    }
    else if (state2 == 0) {
        trans = (double) Fsgn * tmp_trans;
    }
    else return 0;

    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, trans) shared(tmp_v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {
                dmv = trans * tmp_v1buf[j + 1];
                tmp_v0[j + 1] += dmv;
                dam_pr += conj(tmp_v1[j + 1]) * dmv;
        }
    return (dam_pr);
#endif
}/*double complex X_GC_Cis_MPI*/


/**
  *
  * Single creation/annihilation operator
  * in the inter process region for HubbardGC.
  *
  * @author Mitsuaki Kawamura (The University of Tokyo)
  * @author Kazuyoshi Yoshimi (The University of Tokyo)
  * @author Youhei Yamaji (The University of Tokyo)
  */
double complex X_GC_Ajt_MPI(
				       int org_isite,
				       int org_ispin,
				       double complex tmp_trans,
  double complex *tmp_v0 /**< [out] Result v0 += H v1*/,
  double complex *tmp_v1 /**< [in] v0 += H v1*/,
  unsigned long int idim_max,
  double complex *tmp_v1buf,
  long unsigned int *Tpow
  ) {
#ifdef MPI
    int mask2, state2, ierr, origin, bit2diff, Fsgn;
    unsigned long int idim_max_buf, j;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    // org_isite >= Nsite
    mask2 = (int) Tpow[2 * org_isite + org_ispin];

    origin = myrank ^ mask2; // XOR
    state2 = origin & mask2;

    //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
    //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin. 

    bit2diff = myrank - ((2*mask2-1) & myrank);

    //SgnBit((unsigned long int) (origin & bit2diff), &Fsgn); // Fermion sign
    SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign

    ierr = MPI_Sendrecv(&idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        tmp_v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    if (state2 == 0) {
        trans = 0;
    }
    else if (state2 == mask2) {
        trans = (double) Fsgn * tmp_trans;
    }
    else return 0;

    dam_pr = 0.0;
#pragma omp parallel for default(none) reduction(+:dam_pr) private(j, dmv) \
  firstprivate(idim_max_buf, trans) shared(tmp_v1buf, tmp_v1, tmp_v0)
        for (j = 0; j < idim_max_buf; j++) {
                dmv = trans * tmp_v1buf[j + 1];
                tmp_v0[j + 1] += dmv;
                dam_pr += conj(tmp_v1[j + 1]) * dmv;
        }
    return (dam_pr);
#endif
}/*double complex X_GC_Ajt_MPI*/

int X_Cis(
        long unsigned int j,
        long unsigned int is1_spin,
        long unsigned int *tmp_off,
        long unsigned int *list_1_org,
        long unsigned int *list_2_1_target,
        long unsigned int *list_2_2_target,
        const long unsigned int _irght,
        const long unsigned int _ilft,
        const long unsigned int _ihfbit

) {

    long unsigned int list_1_j, list_1_off;
    long unsigned int ibit_tmp_1;
    long unsigned int bit;
    int sgn, ipsgn;

    list_1_j = list_1_org[j];

    ibit_tmp_1 = (list_1_j & is1_spin);
    // is1_spin >= 1
    // is1_spin = Tpow[2*isite + ispin]

    *tmp_off = 0;

    if (ibit_tmp_1 == 0) {
        // able to create an electron at the is1_spin state
        bit = list_1_j - ( list_1_j & (2*is1_spin-1) );
        SgnBit(bit, &sgn); // Fermion sign
        ipsgn = 1;
#ifdef MPI
        SgnBit(myrank, &ipsgn); // Fermion sign
#endif
        list_1_off = list_1_j | is1_spin; // OR

        GetOffComp(list_2_1_target, list_2_2_target, list_1_off, _irght, _ilft, _ihfbit, tmp_off);
        sgn *=ipsgn;
        return (sgn );
    }   else {
            *tmp_off = 1;
            return 0;
        }

}

double complex X_Cis_MPI(
        int org_isite,
        unsigned int org_ispin,
        double complex tmp_trans,
        double complex *tmp_v0,
        double complex *tmp_v1,
        double complex *tmp_v1buf,
        unsigned long int idim_max,
        long unsigned int *Tpow,
        long unsigned int *list_1_org,
        long unsigned int *list_1buf_org,
        long unsigned int *list_2_1_target,
        long unsigned int *list_2_2_target,
        const long unsigned int _irght,
        const long unsigned int _ilft,
        const long unsigned int _ihfbit
) {
#ifdef MPI
    int mask2, state2, ierr, origin, bit2diff, Fsgn;
    unsigned long int idim_max_buf, j, ioff;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    // org_isite >= Nsite
    mask2 = (int) Tpow[2 * org_isite + org_ispin];

    origin = myrank ^ mask2; // XOR
    state2 = origin & mask2;

    //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
    //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

    bit2diff = myrank - ((2*mask2-1) & myrank);

    SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign

    ierr = MPI_Sendrecv(&idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1_org, idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf_org, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        tmp_v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    if (state2 == mask2) {
        trans = 0;
    }
    else if (state2 == 0) {
        trans = (double) Fsgn * tmp_trans;
    }
    else return 0;

    dam_pr = 0.0;
#pragma omp parallel for default(none) private(j, dmv) \
  firstprivate(idim_max_buf, trans, ioff, list_2_1_target, list_2_2_target) shared(tmp_v1buf, tmp_v1, tmp_v0, list_1buf_org, _irght, _ilft, _ihfbit)
    for (j = 1; j <= idim_max_buf; j++) {//idim_max_buf -> original
        GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j],
                   _irght, _ilft, _ihfbit, &ioff);
        dmv = trans * tmp_v1buf[j];
        tmp_v0[ioff] += dmv;
    }
    return (dam_pr);
#endif
}/*double complex X_GC_Cis_MPI*/


double complex X_Ajt_MPI(
        int org_isite,
        unsigned int org_ispin,
        double complex tmp_trans,
        double complex *tmp_v0,
        double complex *tmp_v1,
        double complex *tmp_v1buf,
        unsigned long int idim_max,
        long unsigned int *Tpow,
        long unsigned int *list_1_org,
        long unsigned int *list_1buf_org,
        long unsigned int *list_2_1_target,
        long unsigned int *list_2_2_target,
        const long unsigned int _irght,
        const long unsigned int _ilft,
        const long unsigned int _ihfbit
){
#ifdef MPI
    int mask2, state2, ierr, origin, bit2diff, Fsgn;
    unsigned long int idim_max_buf, j, ioff;
    MPI_Status statusMPI;
    double complex trans, dmv, dam_pr;

    // org_isite >= Nsite
    mask2 = (int) Tpow[2 * org_isite + org_ispin];

    origin = myrank ^ mask2; // XOR
    state2 = origin & mask2;

    //if state2 = mask2, the state (org_isite, org_ispin) is not occupied in myrank
    //origin: if the state (org_isite, org_ispin) is occupied in myrank, the state is not occupied in origin.

    bit2diff = myrank - ((2*mask2-1) & myrank);

    SgnBit((unsigned long int) (bit2diff), &Fsgn); // Fermion sign
    printf("Debug0: rank=%d test\n", myrank);

    ierr = MPI_Sendrecv(&idim_max, 1, MPI_UNSIGNED_LONG, origin, 0,
                        &idim_max_buf, 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(list_1_org, idim_max + 1, MPI_UNSIGNED_LONG, origin, 0,
                        list_1buf_org, idim_max_buf + 1, MPI_UNSIGNED_LONG, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    ierr = MPI_Sendrecv(tmp_v1, idim_max + 1, MPI_DOUBLE_COMPLEX, origin, 0,
                        tmp_v1buf, idim_max_buf + 1, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &statusMPI);
    if (ierr != 0) exitMPI(-1);

    if (state2 == 0) {
        trans = 0;
    }
    else if (state2 == mask2) {
        trans = (double) Fsgn * tmp_trans;
    }
    else return 0;

    dam_pr = 0.0;
#pragma omp parallel for default(none) private(j, dmv) \
  firstprivate(idim_max_buf, trans, ioff, list_2_1_target, list_2_2_target) shared(tmp_v1buf, tmp_v1, tmp_v0, list_1buf_org, _irght, _ilft, _ihfbit)
    for (j = 1; j <= idim_max_buf; j++) {
        GetOffComp(list_2_1_target, list_2_2_target, list_1buf_org[j],
                   _irght, _ilft, _ihfbit, &ioff);
        dmv = trans * tmp_v1buf[j];
        tmp_v0[ioff] += dmv;
    }
    return (dam_pr);
#endif
}

double complex X_Ajt(
        long unsigned int j,
        long unsigned int is1_spin,
        long unsigned int *tmp_off,
        long unsigned int *list_1_org,
        long unsigned int *list_2_1_target,
        long unsigned int *list_2_2_target,
        const long unsigned int _irght,
        const long unsigned int _ilft,
        const long unsigned int _ihfbit
) {
    long unsigned int list_1_j, list_1_off;
    long unsigned int ibit_tmp_1;
    long unsigned int bit;
    int sgn, ipsgn;

    list_1_j = list_1_org[j];

    ibit_tmp_1 = (list_1_j & is1_spin);
// is1_spin >= 1
// is1_spin = Tpow[2*isite + ispin]

    *tmp_off = 0;
    if (ibit_tmp_1 != 0) {
        // able to delete an electron at the is1_spin state
        bit = list_1_j - ( list_1_j & (2*is1_spin-1) );
        SgnBit(bit, &sgn); // Fermion sign
        ipsgn = 1;
#ifdef MPI
        SgnBit(myrank, &ipsgn); // Fermion sign
#endif        
        list_1_off = list_1_j ^ is1_spin;
        GetOffComp(list_2_1_target, list_2_2_target, list_1_off, _irght, _ilft, _ihfbit, tmp_off);
        sgn *= ipsgn;
        return(sgn);
    } else {
        *tmp_off = 1;
        return 0;
    }

}



int GetSingleExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0, /**< [out] Result v0 = H v1*/
  double complex *tmp_v1 /**< [in] v0 = H v1*/
){

  long int idim_max, idim_maxMPI;
  long unsigned int i,j;
  long unsigned int org_isite,ispin,itype;
  long unsigned int is1_spin;
  int isgn=1;
  double complex tmpphi;
  long unsigned int tmp_off=0;
  //idim_max = X->Check.idim_max;
  idim_max = X->Check.idim_maxOrg;
  //tmp_v0
  if(X->Def.NSingleExcitationOperator == 0){
    return TRUE;
  }
  double complex *tmp_v1bufOrg;
  //set size
#ifdef MPI
    idim_maxMPI = MaxMPI_li(X->Check.idim_maxOrg);
    c_malloc1(tmp_v1bufOrg, idim_maxMPI + 1);
#endif // MPI

  switch(X->Def.iCalcModel){
  case HubbardGC:
    // SingleEx  
    //fprintf(stdoutMPI, "SingleOperation in GetSingleExcitedState Re= %lf ; Im= %lf.\n",
    //creal(X->Def.ParaSingleExcitationOperator[0]),cimag(X->Def.ParaSingleExcitationOperator[0]));
    // X->Def.NSingleExcitationOperator 
    // X->Def.SingleExcitationOperator[0][0]
    // X->Def.ParaSingleExcitationOperator[0]
    // clear all elements of tmp_v0 to zero  
    for(i=0;i<X->Def.NSingleExcitationOperator;i++){
      org_isite = X->Def.SingleExcitationOperator[i][0];
      ispin     = X->Def.SingleExcitationOperator[i][1];
      itype     = X->Def.SingleExcitationOperator[i][2];
      tmpphi    = X->Def.ParaSingleExcitationOperator[i];
      if(itype == 1){
        if( org_isite >= X->Def.Nsite){
          X_GC_Cis_MPI(org_isite,ispin,tmpphi,tmp_v0,tmp_v1,idim_max,tmp_v1bufOrg,X->Def.Tpow);
        }
        else{
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X)	\
  firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_off)
          for(j=1;j<=idim_max;j++){
            is1_spin = X->Def.Tpow[2*org_isite+ispin];
            GC_Cis(j,tmp_v0,tmp_v1,is1_spin,tmpphi,&tmp_off);
          }
        }
      }
      else if(itype == 0){
        if( org_isite >= X->Def.Nsite){
          X_GC_Ajt_MPI(org_isite,ispin,tmpphi,tmp_v0,tmp_v1,idim_max,tmp_v1bufOrg,X->Def.Tpow);
        }
        else{
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X)	\
  firstprivate(idim_max, tmpphi, org_isite, ispin) private(j, is1_spin, tmp_off)
          for(j=1;j<=idim_max;j++){
              is1_spin = X->Def.Tpow[2*org_isite+ispin];
              GC_Ajt(j,tmp_v0,tmp_v1,is1_spin,tmpphi,&tmp_off);
          }
        }
      }
    }

    break;
  
  case KondoGC:
  case Hubbard:
  case Kondo:
      idim_max = X->Check.idim_maxOrg;
      for(i=0;i<X->Def.NSingleExcitationOperator;i++){
          org_isite = X->Def.SingleExcitationOperator[i][0];
          ispin     = X->Def.SingleExcitationOperator[i][1];
          itype     = X->Def.SingleExcitationOperator[i][2];
          tmpphi    = X->Def.ParaSingleExcitationOperator[i];
          is1_spin = X->Def.Tpow[2*org_isite+ispin];
          if(itype == 1){
              if( org_isite >= X->Def.Nsite){
                 X_Cis_MPI(org_isite,ispin,tmpphi,tmp_v0,tmp_v1, tmp_v1bufOrg ,idim_max,\
                 X->Def.Tpow, list_1_org, list_1buf_org, list_2_1, list_2_2, \
                 X->Large.irght, X->Large.ilft,X->Large.ihfbit);
              }
              else{
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X, list_1_org)	\
  firstprivate(idim_max, tmpphi, org_isite, ispin, list_2_1, list_2_2, is1_spin) private(j,  isgn,tmp_off)
                  for(j=1;j<=idim_max;j++){//idim_max -> original dimension
                      isgn=X_Cis(j, is1_spin, &tmp_off, list_1_org, list_2_1, list_2_2, X->Large.irght, X->Large.ilft,X->Large.ihfbit);
                      tmp_v0[tmp_off] += tmp_v1[j]*isgn*tmpphi;
                  }
              }
          }
          else if(itype == 0){
              if( org_isite >= X->Def.Nsite){
                  X_Ajt_MPI(org_isite,ispin,tmpphi,tmp_v0,tmp_v1, tmp_v1bufOrg, \
                            idim_max,X->Def.Tpow, list_1_org, list_1buf_org, \
                            list_2_1, list_2_2, X->Large.irght, X->Large.ilft,X->Large.ihfbit);
              }
              else{
#pragma omp parallel for default(none) shared(tmp_v0, tmp_v1, X, list_1_org, list_1) \
  firstprivate(idim_max, tmpphi, org_isite, ispin, list_2_1, list_2_2, is1_spin, myrank) private(j, isgn, tmp_off)
                  for(j=1;j<=idim_max;j++){//idim_max -> original dimension
                      isgn=X_Ajt(j, is1_spin, &tmp_off, list_1_org, list_2_1, list_2_2, X->Large.irght, X->Large.ilft,X->Large.ihfbit);
                      tmp_v0[tmp_off] += tmp_v1[j]*isgn*tmpphi;
                  }
              }
          }
      }
      break;

  case Spin:
  case SpinGC:
    return FALSE;
    break;

  default:
    return FALSE;
  }

  return TRUE;
}

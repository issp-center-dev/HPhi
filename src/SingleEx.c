#include "bitcalc.h"
#include "SingleEx.h"
#include "SingleExHubbard.h"

///
/// Calculation of single excited state
/// Target System: Hubbard, Kondo
/// \param X define list to get and put information of calcuation
/// \param tmp_v0 [out] Result v0 = H v1
/// \param tmp_v1 [in] v0 = H v1
/// \returns TRUE: Normally finished
/// \returns FALSE: Abnormally finished
/// \authour Kazuyoshi Yoshimi
/// \version 1.2
int GetSingleExcitedState
        (
                struct BindStruct *X,
                double complex *tmp_v0,
                double complex *tmp_v1
        ){
    int iret=0;
    //tmp_v0
    if(X->Def.NSingleExcitationOperator == 0){
        return TRUE;
    }

    switch(X->Def.iCalcModel){
        case HubbardGC:
            iret=GetSingleExcitedStateHubbardGC(X,tmp_v0, tmp_v1);
            break;

        case KondoGC:
        case Hubbard:
        case Kondo:
            iret=GetSingleExcitedStateHubbard(X,tmp_v0, tmp_v1);
            break;

        case Spin:
        case SpinGC:
            iret=FALSE;
            break;

        default:
            iret=FALSE;
            break;
    }

    return iret;
}

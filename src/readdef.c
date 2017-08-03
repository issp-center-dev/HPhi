/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*-------------------------------------------------------------
 *[ver.2008.11.4]
 *  Read Definition files
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 *-------------------------------------------------------------*/


/*=================================================================================================*/

/**
 * @file   readdef.c
 * 
 * @brief  File to define functions of reading input files
 * 
 * 
 */


#include "Common.h"
#include "readdef.h"
#include <ctype.h>
#include "LogMessage.h"
#include "wrapperMPI.h"

/**
 * Keyword List in NameListFile.
 **/

static char cKWListOfFileNameList[][D_CharTmpReadDef]={
  "CalcMod",
  "ModPara",
  "LocSpin",
  "Trans",
  "CoulombIntra",
  "CoulombInter",
  "Hund",
  "PairHop",
  "Exchange",
  "InterAll",
  "OneBodyG",
  "TwoBodyG",
  "PairLift",
  "Ising",
  "Boost",
  "SingleExcitation",
  "PairExcitation",
  "SpectrumVec",
  "Laser",
  "TEOneBody",
  "TETwoBody"
};

int D_iKWNumDef = sizeof(cKWListOfFileNameList)/sizeof(cKWListOfFileNameList[0]);

/**
 * File Name List in NameListFile.
 **/
static char (*cFileNameListFile)[D_CharTmpReadDef];

int CheckTETransferHermite(struct DefineList *X,
                           const int NTETrans,
                           const int idx);


int CheckInterAllCondition(
        int iCalcModel,
        int Nsite,
        int iFlgGeneralSpin,
        int *iLocSpin,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4
);

int InputInterAllInfo(
        int *icnt_interall,
        int **iInterAllInfo,
        double complex *cInterAllValue,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4,
        double re_value, double im_value
);


/**
 * @brief Error Function of reading def files.
 * @param[in] _defname name of def file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileError(
                     const	char *defname
                     ){
  fprintf(stdoutMPI, cErrReadDefFile, defname);
  return (-1);
}

/**
 * @brief Function of Validating value.
 * @param[in] icheckValue value to validate.
 * @param[in] ilowestValue lowest value which icheckValue can be set.
 * @param[in] iHighestValue heighest value which icheckValue can be set.
 * @retval 0 value is correct.
 * @retval -1 value is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ValidateValue(
                  const int icheckValue,
                  const int ilowestValue, 
                  const int iHighestValue
                  ){

  if(icheckValue < ilowestValue || icheckValue > iHighestValue){
    return(-1);
  }
  return 0;
}

/**
 * @brief Function of Checking keyword in NameList file.
 * @param[in] _cKW keyword candidate
 * @param[in] _cKWList Reffercnce of keyword List
 * @param[in] iSizeOfKWidx number of keyword
 * @param[out] _iKWidx index of keyword
 * @retval 0 keyword is correct.
 * @retval -1 keyword is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckKW(
            const char* cKW,
            char  cKWList[][D_CharTmpReadDef],
            int iSizeOfKWidx,
            int* iKWidx
            ){
  *iKWidx=-1;
  int itmpKWidx;
  int iret=-1;
  for(itmpKWidx=0; itmpKWidx<iSizeOfKWidx; itmpKWidx++){
    if(strcmp(cKW,"")==0){
      break;
    }
    else if(CheckWords(cKW, cKWList[itmpKWidx])==0){
      iret=0;
      *iKWidx=itmpKWidx;
    }
  }
  return iret;
}

/**
 * @brief Function of Getting keyword and it's variable from characters.
 * @param[in] _ctmpLine characters including keyword and it's variable 
 * @param[out] _ctmp keyword
 * @param[out] _itmp variable for a keyword
 * @retval 0 keyword and it's variable are obtained.
 * @retval 1 ctmpLine is a comment line.
 * @retval -1 format of ctmpLine is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int GetKWWithIdx(
                 char *ctmpLine,
                 char *ctmp,
                 int *itmp
                 )
{
  char *ctmpRead;
  char *cerror;
  char csplit[] = " ,.\t\n";
  if(*ctmpLine=='\n') return 1;
  ctmpRead = strtok(ctmpLine, csplit);
  if(strncmp(ctmpRead, "=", 1)==0 || strncmp(ctmpRead, "#", 1)==0 || ctmpRead==NULL){
    return 1;
  }
  strcpy(ctmp, ctmpRead);
    
  ctmpRead = strtok( NULL, csplit );
  *itmp = strtol(ctmpRead, &cerror, 0);
  //if ctmpRead is not integer type
  if(*cerror != '\0'){
    fprintf(stdoutMPI, cErrDefFileFormat, cerror);
    return(-1);
  }

  ctmpRead = strtok( NULL, csplit );
  if(ctmpRead != NULL){
    fprintf(stdoutMPI, cErrDefFileFormat, ctmpRead);
    return(-1);
  }
    
  return 0;
}

/**
 * @brief Function of Reading calcmod file.
 * @param[in]  _defname file name to read.
 * @param[out] X Define List for getting flags of calc-mode.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ReadcalcmodFile(
                    const char *defname,
                    struct DefineList *X
                    )
{
  FILE *fp;
  int itmp, iret;
  char ctmpLine[D_CharTmpReadDef+D_CharKWDMAX];
  char ctmp[D_CharKWDMAX];
  X->iCalcType=0;
  X->iFlgFiniteTemperature=0;
  X->iCalcModel=0;
  X->iOutputMode=0;
  X->iCalcEigenVec=0;
  X->iInitialVecType=0;
  X->iOutputEigenVec=0;
  X->iInputEigenVec=0;
  X->iOutputHam=0;
  X->iInputHam=0;
  X->iFlgCalcSpec=0;
  X->iReStart=0;
  X->iFlgMPI=0;

  /*=======================================================================*/
  fp = fopenMPI(defname, "r");
  if(fp==NULL) return ReadDefFileError(defname);
  /* read Parameters from calcmod.def*/
  while( fgetsMPI(ctmpLine, D_CharTmpReadDef+D_CharKWDMAX, fp)!=NULL ){
    if( (iret=GetKWWithIdx(ctmpLine, ctmp, &itmp)) !=0){
      if(iret==1) continue;
      return(-1);
    }   
    if(CheckWords(ctmp, "CalcType")==0){
      X->iCalcType=itmp;
    }
    else if(CheckWords(ctmp, "FlgFiniteTemperature")==0){
      X->iFlgFiniteTemperature = itmp;
    }
    else if(CheckWords(ctmp, "CalcModel")==0){
      X->iCalcModel=itmp;
    }
    else if(CheckWords(ctmp, "OutputMode")==0){
      X->iOutputMode=itmp;
    }
    else if(CheckWords(ctmp, "CalcEigenVec")==0){
      X->iCalcEigenVec=itmp;
    }
    else if(CheckWords(ctmp, "InitialVecType")==0){
      X->iInitialVecType=itmp;
    }
    else if(CheckWords(ctmp, "OutputEigenVec")==0 || CheckWords(ctmp, "OEV")==0){
      X->iOutputEigenVec=itmp;
    }
    else if(CheckWords(ctmp, "InputEigenVec")==0 || CheckWords(ctmp, "IEV")==0){
      X->iInputEigenVec=itmp;
    }
    else if(CheckWords(ctmp, "OutputHam")==0){
      X->iOutputHam=itmp;
    }
    else if(CheckWords(ctmp, "InputHam")==0){
      X->iInputHam=itmp;
    }
    else if(CheckWords(ctmp, "CalcSpec")==0 || CheckWords(ctmp, "CalcSpectrum")==0){
      X->iFlgCalcSpec=itmp;
    }
    else if(CheckWords(ctmp, "ReStart")==0){
      X->iReStart=itmp;
    }
    else{
      fprintf(stdoutMPI, cErrDefFileParam, defname, ctmp);
      return(-1);
    }
  }
  fclose(fp);
    
  /* Check values*/
  if(ValidateValue(X->iCalcModel, 0, NUM_CALCMODEL-1)){
    fprintf(stdoutMPI, cErrCalcType, defname);
    return (-1);
  }
  if(ValidateValue(X->iCalcType, 0, NUM_CALCTYPE-1)){
    fprintf(stdoutMPI, cErrCalcType, defname);
    return (-1);
  }
  if(ValidateValue(X->iOutputMode, 0, NUM_OUTPUTMODE-1)){
    fprintf(stdoutMPI, cErrOutputMode, defname);
    return (-1);
  }
  
  if(ValidateValue(X->iCalcEigenVec, -1, NUM_CALCEIGENVEC-1)){
    fprintf(stdoutMPI, cErrCalcEigenVec, defname);
    return (-1);
  }
  
  if(ValidateValue(X->iInitialVecType, 0, NUM_SETINITAILVEC-1)){
    fprintf(stdoutMPI, cErrSetIniVec, defname);
    return (-1);
  }

  if(ValidateValue(X->iOutputHam, 0, NUM_OUTPUTHAM-1)){
    fprintf(stdoutMPI, cErrOutputHam, defname);
    return (-1);
  }
  if(ValidateValue(X->iInputHam, 0, NUM_INPUTHAM-1)){
    fprintf(stdoutMPI, cErrInputHam, defname);
    return (-1);
  }
  if(X->iInputHam == 1 && X->iOutputHam==1){
    fprintf(stdoutMPI, cErrInputOutputHam, defname);
    return (-1);
  }

  if(ValidateValue(X->iReStart, 0, NUM_RESTART-1)){
    fprintf(stdoutMPI, cErrRestart, defname);
    return (-1);
  }

  /* In the case of Full Diagonalization method(iCalcType=2)*/
  if(X->iCalcType==2 && ValidateValue(X->iFlgFiniteTemperature, 0, 1)){
    fprintf(stdoutMPI, cErrFiniteTemp, defname);
    return (-1);
  }

  if(X->iCalcType !=2 && X->iOutputHam ==TRUE) {
    fprintf(stdoutMPI, cErrOutputHamForFullDiag, defname);
    return (-1);
  }

  return 0;
}

/**
 * @brief Function of Fitting FileName
 * @param[in]  _cFileListNameFile file for getting names of input files.
 * @param[out] _cFileNameList arrays for getting names of input files.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int GetFileName(
                const char* cFileListNameFile,
                char cFileNameList[][D_CharTmpReadDef]
                )
{
  FILE *fplist;
  int itmpKWidx=-1;
  char ctmpFileName[D_FileNameMaxReadDef];
  char ctmpKW[D_CharTmpReadDef], ctmp2[256];
  int i;
  for(i=0; i< D_iKWNumDef; i++){
    strcpy(cFileNameList[i],"");
  }

  fplist = fopenMPI(cFileListNameFile, "r");
  if(fplist==NULL) return ReadDefFileError(cFileListNameFile);

  while(fgetsMPI(ctmp2, 256, fplist) != NULL){
    memset(ctmpKW, '\0', strlen(ctmpKW));
    memset(ctmpFileName, '\0', strlen(ctmpFileName));
    sscanf(ctmp2,"%s %s\n", ctmpKW, ctmpFileName);

    if(strncmp(ctmpKW, "#", 1)==0 || *ctmp2=='\n' || (strcmp(ctmpKW, "")&&strcmp(ctmpFileName,""))==0){
      continue;
    }
    else if(strcmp(ctmpKW, "")*strcmp(ctmpFileName, "")==0){
      fprintf(stdoutMPI, cErrKW_InCorPair, cFileListNameFile);
      fclose(fplist);
      return(-1);
    }
    /*!< Check KW */
    if( CheckKW(ctmpKW, cKWListOfFileNameList, D_iKWNumDef, &itmpKWidx)!=0 ){
      fprintf(stdoutMPI, cErrKW, ctmpKW, cFileListNameFile);
      fprintf(stdoutMPI, "%s", cErrKW_ShowList);
      for(i=0; i<D_iKWNumDef;i++){
	    fprintf(stdoutMPI, "%s \n", cKWListOfFileNameList[i]);
      }
      fclose(fplist);
      return(-1);
    }
    /*!< Check cFileNameList to prevent from double registering the file name */
    if(strcmp(cFileNameList[itmpKWidx], "") !=0){
      fprintf(stdoutMPI, cErrKW_Same, cFileListNameFile);
      fclose(fplist);
      return(-1);
    }

    /*!< Copy FileName */
    strcpy(cFileNameList[itmpKWidx], ctmpFileName);
  }
  fclose(fplist);  
  return 0;
}

/** 
 * @brief  Function of reading information about "ModPara" file and total number of parameters from other def files.
 *
 * @param[in] _xNameListFile List of Input File names.
 * @param[out] XX Define List for getting flags of calc-mode.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileNInt(
                    char *xNameListFile, 
                    struct DefineList *X,
                    struct BoostList *xBoost
                    )
{
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_CharTmpReadDef], ctmp2[256];
  int i,itmp;
  unsigned int iline=0;
  X->nvec=0;
  X->iFlgSpecOmegaMax=FALSE;
  X->iFlgSpecOmegaMin=FALSE;
  X->iFlgSpecOmegaOrg=FALSE;
  X->iNOmega=1000;
  X->NCond=0;
  X->iFlgSzConserved=FALSE;
  X->dcOmegaOrg=0;
  int iReadNCond=FALSE;
  xBoost->flgBoost=FALSE;	
  InitializeInteractionNum(X);
  NumAve=1;
  X->Param.ExpecInterval=1;
  cFileNameListFile = malloc(sizeof(char)*D_CharTmpReadDef*D_iKWNumDef);

  fprintf(stdoutMPI, cReadFileNamelist, xNameListFile); 
  if(GetFileName(xNameListFile, cFileNameListFile)!=0){
    return(-1);
  }

  /*=======================================================================*/
  int iKWidx=0;
  //Check the existence of Essensial Files.
  X->READ=0;
  X->WRITE=0;

  for(iKWidx=0; iKWidx< D_iKWNumDef; iKWidx++){ 
    strcpy(defname, cFileNameListFile[iKWidx]);
    if(strcmp(defname,"")==0){
      switch (iKWidx){
      case KWCalcMod:
      case KWModPara:
      case KWLocSpin:
        fprintf(stdoutMPI, cErrMakeDef, cKWListOfFileNameList[iKWidx]);
        return(-1);
      default:
        break;
      }
    }
  } 

  
  for(iKWidx=0; iKWidx< D_iKWNumDef; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);

    if (strcmp(defname, "") == 0) continue;
    if(iKWidx==KWSpectrumVec){
      continue;
    }
    fprintf(stdoutMPI, cReadFile, defname, cKWListOfFileNameList[iKWidx]);
    fp = fopenMPI(defname, "r");
    if (fp == NULL) return ReadDefFileError(defname);
    switch (iKWidx) {
      case KWCalcMod:
        /* Read calcmod.def---------------------------------------*/
        if (ReadcalcmodFile(defname, X) != 0) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
            break;

      case KWModPara:
        /* Read modpara.def---------------------------------------*/
        //TODO: add error procedure here when parameters are not enough.
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &itmp); //2
            fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp); //3
            fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp); //4
            fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp); //5
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %s\n", ctmp, X->CDataFileHead); //6
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %s\n", ctmp, X->CParaFileHead); //7
            fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);   //8
            double dtmp, dtmp2;

            X->read_hacker = 1;
            while (fgetsMPI(ctmp2, 256, fp) != NULL) {
              if (*ctmp2 == '\n') continue;
              sscanf(ctmp2, "%s %lf %lf\n", ctmp, &dtmp, &dtmp2);
              if (CheckWords(ctmp, "Nsite") == 0) {
                X->Nsite = (int) dtmp;
              }
              else if (CheckWords(ctmp, "Nup") == 0) {
                X->Nup = (int) dtmp;
              }
              else if (CheckWords(ctmp, "Ndown") == 0) {
                X->Ndown = (int) dtmp;
                X->Total2Sz = X->Nup - X->Ndown;
              }
              else if (CheckWords(ctmp, "2Sz") == 0) {
                X->Total2Sz = (int) dtmp;
                X->iFlgSzConserved = TRUE;
              }
              else if (CheckWords(ctmp, "Ncond") == 0) {
                if((int) dtmp <0) {
                  fprintf(stdoutMPI, cErrNcond, defname);
                  return (-1);
                }
                X->NCond = (int) dtmp;
                iReadNCond = TRUE;
              }
              else if (CheckWords(ctmp, "Lanczos_max") == 0) {
                X->Lanczos_max = (int) dtmp;
              }
              else if (CheckWords(ctmp, "initial_iv") == 0) {
                X->initial_iv = (int) dtmp;
              }
              else if (CheckWords(ctmp, "nvec") == 0) {
                X->nvec = (int) dtmp;
              }
              else if (CheckWords(ctmp, "exct") == 0) {
                X->k_exct = (int) dtmp;
              }
              else if (CheckWords(ctmp, "LanczosEps") == 0) {
                X->LanczosEps = (int) dtmp;
              }
              else if (CheckWords(ctmp, "LanczosTarget") == 0) {
                X->LanczosTarget = (int) dtmp;
              }
              else if (CheckWords(ctmp, "LargeValue") == 0) {
                LargeValue = dtmp;
              }
              else if (CheckWords(ctmp, "NumAve") == 0) {
                NumAve = (int) dtmp;
              }
              else if(strcmp(ctmp, "TimeSlice")==0){
                X->Param.TimeSlice=dtmp;
              }
              else if(strcmp(ctmp, "ExpandCoef")==0){
                X->Param.ExpandCoef=(int)dtmp;
              }
              else if(strcmp(ctmp, "OutputInterval")==0){
                X->Param.OutputInterval=(int)dtmp;
              }
              else if (CheckWords(ctmp, "ExpecInterval") == 0) {
                X->Param.ExpecInterval = (int) dtmp;
              }
              else if(strcmp(ctmp, "Tinit")==0){
                X->Param.Tinit=dtmp;
              }
              else if (CheckWords(ctmp, "CalcHS") == 0) {
                X->read_hacker = (int) dtmp;
              }
              else if(CheckWords(ctmp, "OmegaMax")==0){
                X->dcOmegaMax=dtmp+dtmp2*I;
                X->iFlgSpecOmegaMax=TRUE;
              }
              else if(CheckWords(ctmp, "OmegaMin")==0){
                X->dcOmegaMin =dtmp+dtmp2*I;
                X->iFlgSpecOmegaMin=TRUE;
              }
              else if(CheckWords(ctmp, "OmegaIm")==0){
                X->dcOmegaOrg +=dtmp*I;
                X->iFlgSpecOmegaOrg=TRUE;
              }
                else if(CheckWords(ctmp, "OmegaOrg")==0){
                X->dcOmegaOrg +=dtmp+dtmp2*I;
                X->iFlgSpecOmegaOrg=TRUE;
              }
              else if(CheckWords(ctmp, "NOmega")==0){
                X->iNOmega=(int)dtmp;
              }
              else if(CheckWords(ctmp, "TargetTPQRand")==0) {
                X->irand=(int)dtmp;
              }
              else {
                return (-1);
              }
            }
            break;

      case KWLocSpin:
        // Read locspn.def
        X->iFlgGeneralSpin = FALSE;
            fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NLocSpn));
            break;
      case KWTrans:
        // Read transfer.def
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NTransfer));
            break;
      case KWCoulombIntra:
        /* Read coulombintra.def----------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NCoulombIntra));
            break;
      case KWCoulombInter:
        /* Read coulombinter.def----------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NCoulombInter));
            break;
      case KWHund:
        /* Read hund.def------------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NHundCoupling));
            break;
      case KWPairHop:
        /* Read pairhop.def---------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NPairHopping));
            X->NPairHopping*=2;
            break;
      case KWExchange:
        /* Read exchange.def--------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NExchangeCoupling));
            break;
      case KWIsing:
        /* Read ising.def--------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NIsingCoupling));
            break;
      case KWPairLift:
        /* Read exchange.def--------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NPairLiftCoupling));
            break;
      case KWInterAll:
        /* Read InterAll.def--------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NInterAll));
            break;
      case KWOneBodyG:
        /* Read cisajs.def----------------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NCisAjt));
            break;
      case KWTwoBodyG:
        /* Read cisajscktaltdc.def--------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%s %d\n", ctmp, &(X->NCisAjtCkuAlvDC));
            break;
      case KWLaser:
        /* Read laser.def--------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
        fgetsMPI(ctmp2, 256, fp);
        sscanf(ctmp2,"%s %d\n", ctmp, &(X->NLaser));
        break;

      case KWTEOneBody:
        if(X->iCalcType != TimeEvolution) break;
        /* Read TEOnebody.def--------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
        fgetsMPI(ctmp2, 256, fp);
        sscanf(ctmp2,"%s %d\n", ctmp, &(X->NTETimeSteps));
        fgetsMPI(ctmp2, 256, fp);
        fgetsMPI(ctmp2, 256, fp);
        fgetsMPI(ctmp2, 256, fp);
        int iTETransMax=0;
        if(X->NTETimeSteps>0) {
          while (fgetsMPI(ctmp2, 256, fp) != NULL) {
            sscanf(ctmp2, "%lf %d \n", &dtmp, &itmp);
            for (i = 0; i < itmp; ++i) {
              fgetsMPI(ctmp2, 256, fp);
            }
            if(iTETransMax < itmp) iTETransMax=itmp;
          }
        }
      X->NTETransferMax=iTETransMax;
      break;

      case KWTETwoBody:
        if(X->iCalcType != TimeEvolution) break;
        /* Read TETwobody.def--------------------------------*/
        fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
        fgetsMPI(ctmp2, 256, fp);
        sscanf(ctmp2,"%s %d\n", ctmp, &(X->NTETimeSteps));
        fgetsMPI(ctmp2, 256, fp);
        fgetsMPI(ctmp2, 256, fp);
        fgetsMPI(ctmp2, 256, fp);
        int iTEInterAllMax=0;
        if(X->NTETimeSteps>0) {
          while (fgetsMPI(ctmp2, 256, fp) != NULL) {
            sscanf(ctmp2, "%lf %d \n", &dtmp, &itmp);
            for (i = 0; i < itmp; ++i) {
              fgetsMPI(ctmp2, 256, fp);
            }
            if(iTEInterAllMax < itmp) iTEInterAllMax=itmp;
          }
        }
        X->NTEInterAllMax=iTEInterAllMax;
        break;


      case KWBoost:
        /* Read boost.def--------------------------------*/
        xBoost->NumarrayJ = 0;
            xBoost->W0 = 0;
            xBoost->R0 = 0;
            xBoost->num_pivot = 0;
            xBoost->ishift_nspin = 0;
            xBoost->flgBoost = TRUE;
            //first line is skipped
            fgetsMPI(ctmp2, 256, fp);
            //read numarrayJ
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%d\n", &(xBoost->NumarrayJ));
            //skipp arrayJ
            for (iline = 0; iline < xBoost->NumarrayJ * 3; iline++) {
              fgetsMPI(ctmp2, 256, fp);
            }
            //read W0 R0 num_pivot ishift_nspin
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%ld %ld %ld %ld\n", &(xBoost->W0), &(xBoost->R0), &(xBoost->num_pivot),
                   &(xBoost->ishift_nspin));

            break;

    case KWSingleExcitation:
      /* Read singleexcitation.def----------------------------------------*/
      fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fgetsMPI(ctmp2, 256, fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &(X->NSingleExcitationOperator));
      break;

    case KWPairExcitation:
      /* Read pairexcitation.def----------------------------------------*/
      fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fgetsMPI(ctmp2, 256, fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &(X->NPairExcitationOperator));
      break;

    default:
      fprintf(stdoutMPI, "%s", cErrIncorrectDef);
      fclose(fp);
      return (-1);
    }
    /*=======================================================================*/
    fclose(fp);
  }

  //Sz, Ncond
  switch(X->iCalcModel){
  case Spin:
  case Hubbard:
  case Kondo: 
  case SpinlessFermion:
   
    if(iReadNCond==TRUE){
      if(X->iCalcModel==Spin){
        fprintf(stdoutMPI, "For Spin, Ncond should not be defined.\n");
        return(-1);
      }
      else{
        if(X->iFlgSzConserved==TRUE){
          if(X->iCalcModel==SpinlessFermion){
            fprintf(stdoutMPI, "  Warning: For Spinless fermion, 2Sz should not be defined.\n");
            X->Ne=X->NCond;  
            X->Nup=X->NCond;
            X->Ndown=0;
            break;
          }
          X->Nup=X->NLocSpn+X->NCond+X->Total2Sz;
          X->Ndown=X->NLocSpn+X->NCond-X->Total2Sz;
          X->Nup/=2;
          X->Ndown/=2;
        }
        else{
          if(X->iCalcModel == Hubbard){
            X->Ne=X->NCond;
            if(X->Ne <1){
              fprintf(stdoutMPI, "Ncond is incorrect.\n");
              return(-1);
            }
            X->iCalcModel=HubbardNConserved;
          }
          else if(X->iCalcModel ==SpinlessFermion){
            X->Ne=X->NCond;  
            X->Nup=X->NCond;
            X->Ndown=0;
          }
          else{
            fprintf(stdoutMPI, " 2Sz is not defined.\n");
            return(-1);
          }
        }
      }
    }
    else if(iReadNCond == FALSE && X->iFlgSzConserved==TRUE){
      if(X->iCalcModel != Spin){
        fprintf(stdoutMPI, " NCond is not defined.\n");
        return(-1);
      }
      X->Nup=X->NLocSpn+X->Total2Sz;
      X->Ndown=X->NLocSpn-X->Total2Sz;
      X->Nup /= 2;
      X->Ndown /= 2;
    }
    else{
      if(X->Nup==0 && X->Ndown==0){
        if(X->iCalcModel == Spin){
          fprintf(stdoutMPI, " 2Sz is not defined.\n");
          return(-1);
        }
        else{
          fprintf(stdoutMPI, " NCond is not defined.\n");
          return(-1);
        }
      }
    }
    
    if(X->iCalcModel == Spin){
      X->Ne=X->Nup;
    }
    else{
      if(X->Ne==0) {
        X->Ne = X->Nup + X->Ndown;
      }
      if(X->NLocSpn>X->Ne){
        fprintf(stdoutMPI, "%s", cErrNLoc);
        fprintf(stdoutMPI, "NLocalSpin=%d, Ne=%d\n", X->NLocSpn, X->Ne);
        return(-1);
      }
    }
    break;
  case SpinGC:
  case KondoGC:
  case HubbardGC:
  case SpinlessFermionGC:  
    if(iReadNCond == TRUE || X->iFlgSzConserved ==TRUE){
      fprintf(stdoutMPI, "\n  Warning: For GC, both Ncond and 2Sz should not be defined.\n");
      //return(-1);
    }
    break;
  default:
    break;
  }

    /* Check values (Positive)*/
    if(X->Nsite<=0) {// Nsite must be positve
      fprintf(stdoutMPI, cErrNsite, defname);
      return (-1);
    }
    if(X->Lanczos_max<=0) {// Lanczos_max must be positive
      fprintf(stdoutMPI, cErrLanczos_max, defname);
      return (-1);
    }
    if(X->LanczosEps<=0) {// Lanczos_eps must be positive
      fprintf(stdoutMPI, cErrLanczos_eps, defname);
      return (-1);
    }
    if(NumAve<=0) { // Average number must be positive
      fprintf(stdoutMPI, cErrNumAve, defname);
      return (-1);
    }
    if(X->Param.ExpecInterval<=0){// Interval to calculate expected values must be positive
      fprintf(stdoutMPI, cErrExpecInterval, defname);
      return (-1);
    }
    if(X->nvec==0){
      X->nvec=X->Lanczos_max;
    }

    if(X->nvec < X->k_exct){
        X->nvec=X->k_exct;
    }
    if (X->LanczosTarget < X->k_exct) X->LanczosTarget = X->k_exct;

    if(ValidateValue(X->k_exct, 1, X->nvec)) {
      fprintf(stdoutMPI, cErrLanczosExct, defname, X->nvec);
      return (-1);
    }

    if(ValidateValue(X->LanczosTarget, X->k_exct, X->Lanczos_max)){
      fprintf(stdoutMPI, cErrLanczosTarget, defname, X->LanczosTarget, X->k_exct, X->Lanczos_max);
      return (-1);
    }

  X->fidx = 0;
  X->NeMPI=X->Ne;
  X->NupMPI=X->Nup;
  X->NdownMPI=X->Ndown;
  X->NupOrg=X->Nup;
  X->NdownOrg=X->Ndown;
  return 0;
}

/** 
 * @brief function of reading def files to get keyword index
 * 
 * @param X define list to get and put informations for calcuation
 * 
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileIdxPara(
                       struct DefineList *X,
                       struct BoostList *xBoost
                       )
{
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_CharTmpReadDef], ctmp2[256];

  unsigned int i,idx, itype;
  int xitmp[8];
  int iKWidx=0;
  int iboolLoc=0;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  double dvalue_re, dvalue_im;
  double dArrayValue_re[3]; 
  int icnt_diagonal=0;
  int ieps_CheckImag0=-12;
  eps_CheckImag0=pow(10.0, ieps_CheckImag0);
  unsigned int iline=0;
  int ilineIn=0;
  int ilineIn2=0;
  int itmp=0;
  int icnt_trans=0;
  int iflg_trans=0;
  int icnt_interall=0;
  int iflg_interall=0;

  unsigned int iloop=0;

  for(iKWidx=KWLocSpin; iKWidx< D_iKWNumDef; iKWidx++){
    strcpy(defname, cFileNameListFile[iKWidx]);
    if(strcmp(defname,"")==0 || iKWidx==KWSpectrumVec) continue;
    fprintf(stdoutMPI, cReadFileNamelist, defname);
    fp = fopenMPI(defname, "r");
    if(fp==NULL) return ReadDefFileError(defname);
    if(iKWidx != KWBoost){
      for(i=0;i<IgnoreLinesInDef;i++) fgetsMPI(ctmp, sizeof(ctmp)/sizeof(char), fp);
    }
    
    idx=0;    
    /*=======================================================================*/
    switch(iKWidx){
    case KWLocSpin:
      /* Read locspn.def----------------------------------------*/
      while( fgetsMPI(ctmp2, 256, fp) != NULL){
        if(idx==X->Nsite){
          fclose(fp);
          return ReadDefFileError(defname);
        }

        sscanf(ctmp2, "%d %d\n", &(xitmp[0]), &(xitmp[1]) );
        X->LocSpn[xitmp[0]] = xitmp[1];
        X->SiteToBit[xitmp[0]]=(X->LocSpn[xitmp[0]]+1);//2S+1
        if(CheckSite(xitmp[0], X->Nsite) !=0){
          fclose(fp);
          return ReadDefFileError(defname);
        }       
        idx++;
      }
      if(CheckLocSpin(X)==FALSE){
        fclose(fp);
        return ReadDefFileError(defname);
      }

      break;
      
    case KWTrans:
      /* transfer.def--------------------------------------*/
      if(X->NTransfer>0){
        icnt_trans=0;
        while( fgetsMPI(ctmp2, 256, fp) != NULL )
          {
            if(idx==X->NTransfer){
              fclose(fp);
              return ReadDefFileError(defname);
            }

            sscanf(ctmp2, "%d %d %d %d %lf %lf\n",
                   &isite1,
                   &isigma1,
                   &isite2,
                   &isigma2,
                   &dvalue_re,
                   &dvalue_im
                   );

            if(CheckPairSite(isite1, isite2,X->Nsite) !=0){
              fclose(fp);
              return ReadDefFileError(defname);
            }
            
            if(isite1==isite2 && isigma1==isigma2){
              if(fabs(dvalue_im)> eps_CheckImag0){
                //NonHermite
                fprintf(stdoutMPI, cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, dvalue_re, dvalue_im);
                fclose(fp);
                return ReadDefFileError(defname);
              }
            }

            if(X->iCalcModel==Spin){
              if(isite1 != isite2){
                iboolLoc=1;
                fprintf(stdoutMPI, cWarningIncorrectFormatForSpin2, isite1, isite2);
              }
            }	    
            else if(X->iCalcModel==Kondo){
              if(X->LocSpn[isite1]!=ITINERANT || X->LocSpn[isite2] !=ITINERANT){
                if(isite1 != isite2){
                  iboolLoc=1;
                  fprintf(stdoutMPI, cErrIncorrectFormatForKondoTrans, isite1, isite2);
                }
              }
            }
            else if(X->iCalcModel==SpinlessFermion || X->iCalcModel==SpinlessFermionGC){
              if(isigma1 != 0 || isigma2 !=0){
                //Not allowed
                fprintf(stderr, cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, dvalue_re, dvalue_im);
                fclose(fp);
                return ReadDefFileError(defname);
              }
            }
            
            iflg_trans=0;
            for( i=0; i < icnt_trans; i++){
              if(isite1 ==X->GeneralTransfer[i][0] && isite2 == X->GeneralTransfer[i][2]
                 && isigma1 == X->GeneralTransfer[i][1] && isigma2 == X->GeneralTransfer[i][3])
                {
                  X->ParaGeneralTransfer[i] += dvalue_re+dvalue_im*I;
                  iflg_trans=1;
                  continue;
                }
            }
            
            if(iflg_trans == 0){
              X->GeneralTransfer[icnt_trans][0]=isite1;
              X->GeneralTransfer[icnt_trans][1]=isigma1;
              X->GeneralTransfer[icnt_trans][2]=isite2;
              X->GeneralTransfer[icnt_trans][3]=isigma2;
              X->ParaGeneralTransfer[icnt_trans] = dvalue_re+dvalue_im*I;
              icnt_trans++;
            }
            idx++;
          }
	
        if(iboolLoc ==1){
          fclose(fp);
          return(-1);
        }
      }
      
      X->NTransfer = icnt_trans;
      
      if(CheckSpinIndexForTrans(X)==FALSE){
        fclose(fp);
        return(-1);
      }
      
      if(CheckTransferHermite(X) !=0){
        fprintf(stdoutMPI, "%s", cErrNonHermiteTransForAll);
        fclose(fp);
        return(-1);
      }
      break;
      
    case KWCoulombIntra:
      /*coulombintra.def----------------------------------*/
      if(X->NCoulombIntra>0){
        while(fgetsMPI(ctmp2, 256, fp) != NULL){
          if(idx==X->NCoulombIntra){
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %lf\n",
                 &(X->CoulombIntra[idx][0]),
                 &(X->ParaCoulombIntra[idx])
                 );
	  
          if(CheckSite(X->CoulombIntra[idx][0], X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }
          idx++;
        }
      }
      break;

    case KWCoulombInter:
      /*coulombinter.def----------------------------------*/
      if(X->NCoulombInter>0){
        while(fgetsMPI(ctmp2, 256, fp) != NULL){
          if(idx==X->NCoulombInter){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
                 &(X->CoulombInter[idx][0]),
                 &(X->CoulombInter[idx][1]),
                 &(X->ParaCoulombInter[idx])
                 );

          if(CheckPairSite(X->CoulombInter[idx][0], X->CoulombInter[idx][1],X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;

    case KWHund:
      /*hund.def------------------------------------------*/
      if(X->NHundCoupling>0){
        while(fgetsMPI(ctmp2,256,fp) != NULL)
          {
            if(idx==X->NHundCoupling){
              fclose(fp);
              return ReadDefFileError(defname);
            }

            sscanf(ctmp2, "%d %d %lf\n",
                   &(X->HundCoupling[idx][0]),
                   &(X->HundCoupling[idx][1]),
                   &(X->ParaHundCoupling[idx])
                   );

            if(CheckPairSite(X->HundCoupling[idx][0], X->HundCoupling[idx][1],X->Nsite) !=0){
              fclose(fp);
              return ReadDefFileError(defname);
            }

            idx++;
          }
      }
      break;
    case KWPairHop:
      /*pairhop.def---------------------------------------*/
      if(X->iCalcModel == Spin || X->iCalcModel == SpinGC){
        fprintf(stdoutMPI, "PairHop is not active in Spin and SpinGC.\n");
        return(-1);
      }
      
      if(X->NPairHopping>0){
        while(fgetsMPI(ctmp2, 256, fp) != NULL){
          if(idx==X->NPairHopping/2){
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %d %lf\n",
                 &(X->PairHopping[2*idx][0]),
                 &(X->PairHopping[2*idx][1]),
                 &(X->ParaPairHopping[2*idx])
                 );
	  
          if(CheckPairSite(X->PairHopping[2*idx][0], X->PairHopping[2*idx][1],X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }
          X->PairHopping[2*idx+1][0]=X->PairHopping[2*idx][1];
          X->PairHopping[2*idx+1][1]=X->PairHopping[2*idx][0];
          X->ParaPairHopping[2*idx+1]=X->ParaPairHopping[2*idx];
          idx++;
        }
      }
      break;

    case KWExchange:
      /*exchange.def--------------------------------------*/
      if(X->NExchangeCoupling>0){
        while(fgetsMPI(ctmp2,256,fp) != NULL){
          if(idx==X->NExchangeCoupling){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
                 &(X->ExchangeCoupling[idx][0]),
                 &(X->ExchangeCoupling[idx][1]),
                 &(X->ParaExchangeCoupling[idx])
                 );

          if(CheckPairSite(X->ExchangeCoupling[idx][0], X->ExchangeCoupling[idx][1],X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;

    case KWIsing:
      /*ising.def--------------------------------------*/
      if(X->NIsingCoupling>0){
        while(fgetsMPI(ctmp2,256,fp) != NULL){
          if(idx==X->NIsingCoupling){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %lf\n",
                 &isite1,
                 &isite2,
                 &dvalue_re
                 );

          if(CheckPairSite(isite1,isite2,X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          //input into exchange couplings
          X->HundCoupling[X->NHundCoupling+idx][0]=isite1;
          X->HundCoupling[X->NHundCoupling+idx][1]=isite2;
          X->ParaHundCoupling[X->NHundCoupling+idx]= -dvalue_re/2.0;
          //input into inter Coulomb
          X->CoulombInter[X->NCoulombInter+idx][0]=isite1;
          X->CoulombInter[X->NCoulombInter+idx][1]=isite2;
          X->ParaCoulombInter[X->NCoulombInter+idx]=-dvalue_re/4.0;
          idx++;
        }
      }
      break;
      
    case KWPairLift:
      /*pairlift.def--------------------------------------*/
      if(X->NPairLiftCoupling>0){
        if(X->iCalcModel != SpinGC){
          fprintf(stdoutMPI, "PairLift is active only in SpinGC.\n");
          return(-1);
        }
        while(fgetsMPI(ctmp2,256,fp) != NULL)
          {
            if(idx==X->NPairLiftCoupling){
              fclose(fp);
              return ReadDefFileError(defname);
            }

            sscanf(ctmp2, "%d %d %lf\n",
                   &(X->PairLiftCoupling[idx][0]),
                   &(X->PairLiftCoupling[idx][1]),
                   &(X->ParaPairLiftCoupling[idx])
                   );

            if(CheckPairSite(X->PairLiftCoupling[idx][0], X->PairLiftCoupling[idx][1],X->Nsite) !=0){
              fclose(fp);
              return ReadDefFileError(defname);
            }

            idx++;
          }
      }
      break;
      
    case KWInterAll:
      /*interall.def---------------------------------------*/
      X->NInterAll_Diagonal=0;
      X->NInterAll_OffDiagonal=0;
      if(X->NInterAll>0) {
        icnt_interall = 0;
        while (fgetsMPI(ctmp2, 256, fp) != NULL) {
          if (idx == X->NInterAll) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %d %d %d %d %d %d %d %lf %lf\n",
                 &isite1,
                 &isigma1,
                 &isite2,
                 &isigma2,
                 &isite3,
                 &isigma3,
                 &isite4,
                 &isigma4,
                 &dvalue_re,
                 &dvalue_im
          );

          if (CheckInterAllCondition(X->iCalcModel, X->Nsite, X->iFlgGeneralSpin, X->LocSpn,
                                     isite1, isigma1, isite2, isigma2,
                                     isite3, isigma3, isite4, isigma4) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          if (InputInterAllInfo(&icnt_interall,
                                X->InterAll,
                                X->ParaInterAll,
                                isite1, isigma1,
                                isite2, isigma2,
                                isite3, isigma3,
                                isite4, isigma4,
                                dvalue_re, dvalue_im
          ) != 0) {
            icnt_diagonal += 1;
          }
          idx++;
        }
      }

      X->NInterAll = icnt_interall;
      X->NInterAll_Diagonal=icnt_diagonal;
      X->NInterAll_OffDiagonal = X->NInterAll-X->NInterAll_Diagonal;
      /*
        int **InterAll,
        complex double *ParaInterAll,
        const int NInterAll,
        int **InterAllDiagonal,
        double *ParaInterAllDiagonal,
        int **InterAllOffDiagonal,
        complex double *ParaInterAllOffDiagonal,
        int *Chemi,
        int *SpinChemi,
        double *ParaChemi,
        int *NChemi,
        const int iCalcModel
        */
        /*
      if(GetDiagonalInterAll(X)!=0){
        fclose(fp);
        return(-1);
      }
      */

        if(GetDiagonalInterAll(
                X->InterAll, X->ParaInterAll, X->NInterAll,
                X->InterAll_Diagonal, X->ParaInterAll_Diagonal,
                X->InterAll_OffDiagonal, X->ParaInterAll_OffDiagonal,
                X->EDChemi, X->EDSpinChemi, X->EDParaChemi, &X->EDNChemi,
                X->iCalcModel
        )!=0){
          fclose(fp);
          return(-1);
        }

        if(CheckInterAllHermite(X)!=0){
        fprintf(stdoutMPI, "%s", cErrNonHermiteInterAllForAll);
        fclose(fp);
        return(-1);
      }
      break;
      
    case KWOneBodyG:
      /*cisajs.def----------------------------------------*/
      if(X->NCisAjt>0){
        while(fgetsMPI(ctmp2, 256, fp) != NULL){
          if(idx==X->NCisAjt){
            fclose(fp);
            return ReadDefFileError(defname);
          }
          sscanf(ctmp2, "%d %d %d %d\n",
                 &isite1,
                 &isigma1,
                 &isite2,
                 &isigma2);

          if(X->iCalcModel == Spin){
            if(isite1 != isite2){
              fprintf(stdoutMPI, cWarningIncorrectFormatForSpin2, isite1, isite2);
              X->NCisAjt--;
              continue;
            }
          }

          X->CisAjt[ idx ][0] = isite1;
          X->CisAjt[ idx ][1] = isigma1;
          X->CisAjt[ idx ][2] = isite2;
          X->CisAjt[ idx ][3] = isigma2;	  

          if(CheckPairSite(isite1, isite2,X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          idx++;
        }
      }
      break;
      
    case KWTwoBodyG:
      /*cisajscktaltdc.def--------------------------------*/
      if(X->NCisAjtCkuAlvDC>0){
        while(fgetsMPI(ctmp2, 256, fp) != NULL){
          if(idx==X->NCisAjtCkuAlvDC){
            fclose(fp);
            return ReadDefFileError(defname);
          }

          sscanf(ctmp2, "%d %d %d %d %d %d %d %d\n",
                 &isite1,
                 &isigma1,
                 &isite2,
                 &isigma2,
                 &isite3,
                 &isigma3,
                 &isite4,
                 &isigma4
                 );

          if(X->iCalcModel == Spin || X->iCalcModel == SpinGC){
            if(!CheckFormatForSpinInt(isite1, isite2, isite3, isite4)==0){
              X->NCisAjtCkuAlvDC--;
              continue;
            }
          }


          X->CisAjtCkuAlvDC[idx][0] = isite1;
          X->CisAjtCkuAlvDC[idx][1] = isigma1;
          X->CisAjtCkuAlvDC[idx][2] = isite2;
          X->CisAjtCkuAlvDC[idx][3] = isigma2;
          X->CisAjtCkuAlvDC[idx][4] = isite3;
          X->CisAjtCkuAlvDC[idx][5] = isigma3;
          X->CisAjtCkuAlvDC[idx][6] = isite4;
          X->CisAjtCkuAlvDC[idx][7] = isigma4;

          if(CheckQuadSite(isite1, isite2, isite3, isite4,X->Nsite) !=0){
            fclose(fp);
            return ReadDefFileError(defname);
          }
          idx++;
        }
      }
      break;

      case KWLaser:
        //printf("KWLaser\n");
        /*laser.def----------------------------------*/
        if(X->NLaser>0){
          //printf("Read Start\n");
          while(fgetsMPI(ctmp2, 256, fp) != NULL){
            sscanf(ctmp2, "%s %lf\n", &(ctmp), &(X->ParaLaser[idx]));
            //printf("[%d]:%f\n",idx,X->ParaLaser[idx]);
            idx++;
          }
          if(idx!=X->NLaser){
            fclose(fp);
            return ReadDefFileError(defname);
          }
        }
        break;

      case KWTEOneBody:
        if(X->NTETimeSteps>0){
          idx=0;
          while(fgetsMPI(ctmp2, 256, fp) != NULL){
            sscanf(ctmp2, "%lf %d\n", &(X->TETime[idx]), &(X->NTETransfer[idx]));
            for(i=0; i<X->NTETransfer[idx]; ++i ){
              fgetsMPI(ctmp2, 256, fp);
              sscanf(ctmp2, "%d %d %d %d %lf %lf\n",
                     &isite1,
                     &isigma1,
                     &isite2,
                     &isigma2,
                     &dvalue_re,
                     &dvalue_im
                    );
              X->TETransfer[idx][i][0]= isite1;
              X->TETransfer[idx][i][1]= isigma1;
              X->TETransfer[idx][i][2]= isite2;
              X->TETransfer[idx][i][3] = isigma2;
              X->ParaTETransfer[idx][i]=dvalue_re+dvalue_im*I;
            }
            //check Transfer Hermite
            if(CheckTETransferHermite(X, X->NTETransfer[idx], idx)!=0){
              fclose(fp);
              return ReadDefFileError(defname);
            }
            idx++;
          }
          if(idx!=X->NTETimeSteps){
            fclose(fp);
            return ReadDefFileError(defname);
          }
        }
        break;

      case KWTETwoBody:
        if(X->NTETimeSteps>0){
          idx=0;
          while(fgetsMPI(ctmp2, 256, fp) != NULL){
            sscanf(ctmp2, "%lf %d\n", &(X->TETime[idx]), &(X->NTEInterAll[idx]));
            for(i=0; i<X->NTEInterAll[idx]; ++i ){
              fgetsMPI(ctmp2, 256, fp);
              sscanf(ctmp2, "%d %d %d %d %d %d %d %d %lf %lf\n",
                     &isite1,
                     &isigma1,
                     &isite2,
                     &isigma2,
                     &isite3,
                     &isigma3,
                     &isite4,
                     &isigma4,
                     &dvalue_re,
                     &dvalue_im
              );

              if(CheckInterAllCondition(X->iCalcModel, X->Nsite, X->iFlgGeneralSpin,X->LocSpn,
                                        isite1, isigma1, isite2, isigma2,
                                        isite3, isigma3, isite4, isigma4) !=0){
                fclose(fp);
                return ReadDefFileError(defname);
              }

              if (InputInterAllInfo(&icnt_interall,
                                    X->TEInterAll[idx],
                                    X->ParaTEInterAll[idx],
                                    isite1, isigma1,
                                    isite2, isigma2,
                                    isite3, isigma3,
                                    isite4, isigma4,
                                    dvalue_re, dvalue_im
              ) != 0) {
                icnt_diagonal += 1;
              }
            }

            X->NTEInterAll[idx] = icnt_interall;
            X->NTEInterAllDiagonal[idx]=icnt_diagonal;
            X->NTEInterAllOffDiagonal[idx] = icnt_interall-icnt_diagonal;


            if(GetDiagonalInterAll(
                    X->TEInterAll[idx], X->ParaTEInterAll[idx], X->NTEInterAll[idx],
                    X->TEInterAllDiagonal[idx], X->ParaTEInterAllDiagonal[idx],
                    X->TEInterAllOffDiagonal[idx], X->ParaTEInterAllOffDiagonal[idx],
                    X->TEChemi[idx], X->SpinTEChemi[idx], X->ParaTEChemi[idx], &X->NTEChemi[idx],
                    X->iCalcModel
            )!=0){
              fclose(fp);
              return(-1);
            }


            /*
            if(GetTEDiagonalInterAll(X)!=0){
              fclose(fp);
              return(-1);
            }
*/
            if(CheckInterAllHermite(X)!=0){
              fprintf(stdoutMPI, "%s", cErrNonHermiteInterAllForAll);
              fclose(fp);
              return(-1);
            }

            //check Transfer Hermite
            /*
            if(CheckTETransferHermite(X, X->NTETransfer[idx], idx)!=0){
              fclose(fp);
              return ReadDefFileError(defname);
            }
             */
            idx++;
          }

          if(idx!=X->NTETimeSteps){
            fclose(fp);
            return ReadDefFileError(defname);
          }
        }
        break;

      case KWBoost:
      /* boost.def--------------------------------*/
      //input magnetic field
      fgetsMPI(ctmp2, 256, fp);
      sscanf(ctmp2, "%lf %lf %lf\n",
             &dArrayValue_re[0],
             &dArrayValue_re[1],
             &dArrayValue_re[2]);
      for(iline=0; iline<3; iline++){
        xBoost->vecB[iline]= dArrayValue_re[iline];
      }
      
      //this line is skipped;
      fgetsMPI(ctmp2, 256, fp);

      //input arrayJ
      if(xBoost->NumarrayJ>0){
        for(iline=0; iline<xBoost->NumarrayJ; iline++){
          for(ilineIn=0; ilineIn<3; ilineIn++){
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%lf %lf %lf\n",
                   &dArrayValue_re[0],
                   &dArrayValue_re[1],
                   &dArrayValue_re[2]);
            for(ilineIn2=0; ilineIn2<3; ilineIn2++){
              xBoost->arrayJ[iline][ilineIn][ilineIn2]= dArrayValue_re[ilineIn2];
            }	    
          }
        }
      }

      //this line is skipped;
      fgetsMPI(ctmp2, 256, fp);

      //read list_6spin_star
      if(xBoost->num_pivot>0){
        for(iline=0; iline<xBoost->num_pivot; iline++){
          //input
          fgetsMPI(ctmp2, 256, fp);
          sscanf(ctmp2, "%d %d %d %d %d %d %d\n",
                 &xBoost->list_6spin_star[iline][0],
                 &xBoost->list_6spin_star[iline][1],
                 &xBoost->list_6spin_star[iline][2],
                 &xBoost->list_6spin_star[iline][3],
                 &xBoost->list_6spin_star[iline][4],
                 &xBoost->list_6spin_star[iline][5],
                 &xBoost->list_6spin_star[iline][6]
                 ); 
          //copy
          for(iloop=0; iloop<xBoost->R0; iloop++){
            for(itmp=0; itmp<7; itmp++){
              xBoost->list_6spin_star[iloop*xBoost->num_pivot+iline][itmp]=xBoost->list_6spin_star[iline][itmp];
            }
          }   
        }	
      }

      //read list_6spin_pair
      if(xBoost->num_pivot>0){
        for(iline=0; iline<xBoost->num_pivot; iline++){
          //input
          for(ilineIn2=0; ilineIn2<xBoost->list_6spin_star[iline][0]; ilineIn2++){
            fgetsMPI(ctmp2, 256, fp);
            sscanf(ctmp2, "%d %d %d %d %d %d %d\n",
                   &xBoost->list_6spin_pair[iline][0][ilineIn2],
                   &xBoost->list_6spin_pair[iline][1][ilineIn2],
                   &xBoost->list_6spin_pair[iline][2][ilineIn2],
                   &xBoost->list_6spin_pair[iline][3][ilineIn2],
                   &xBoost->list_6spin_pair[iline][4][ilineIn2],
                   &xBoost->list_6spin_pair[iline][5][ilineIn2],
                   &xBoost->list_6spin_pair[iline][6][ilineIn2]
                   ); 

            //copy
            for(iloop=0; iloop<xBoost->R0; iloop++){
              for(itmp=0; itmp<7; itmp++){
                xBoost->list_6spin_pair[iloop*xBoost->num_pivot+iline][itmp][ilineIn2]=xBoost->list_6spin_pair[iline][itmp][ilineIn2];
              }
            }
          }
        }
	
      }

      break;

    case KWSingleExcitation:
      /*singleexcitation.def----------------------------------------*/
      if(X->NSingleExcitationOperator>0) {
        if(X->iCalcModel == Spin || X->iCalcModel == SpinGC) {
          fprintf(stderr, "SingleExcitation is not allowed for spin system.\n");
          fclose(fp);
          return ReadDefFileError(defname);
        }
        while (fgetsMPI(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%d %d %d %lf %lf\n",
                 &isite1,
                 &isigma1,
                 &itype,
                 &dvalue_re,
                 &dvalue_im
                 );

          if (CheckSite(isite1, X->Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }

          X->SingleExcitationOperator[idx][0] = isite1;
          X->SingleExcitationOperator[idx][1] = isigma1;
          X->SingleExcitationOperator[idx][2] = itype;
          X->ParaSingleExcitationOperator[idx] = dvalue_re + I * dvalue_im;
          idx++;
        }
        if (idx != X->NSingleExcitationOperator) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    case KWPairExcitation:
      /*pairexcitation.def----------------------------------------*/
      if(X->NPairExcitationOperator>0) {
        while (fgetsMPI(ctmp2, 256, fp) != NULL) {
          sscanf(ctmp2, "%d %d %d %d %d %lf %lf\n",
                 &isite1,
                 &isigma1,
                 &isite2,
                 &isigma2,
                 &itype,
                 &dvalue_re,
                 &dvalue_im
                 );
          if (CheckPairSite(isite1, isite2, X->Nsite) != 0) {
            fclose(fp);
            return ReadDefFileError(defname);
          }
          if(itype==0){
            X->PairExcitationOperator[idx][0] = isite1;
            X->PairExcitationOperator[idx][1] = isigma1;
            X->PairExcitationOperator[idx][2] = isite2;
            X->PairExcitationOperator[idx][3] = isigma2;
            X->PairExcitationOperator[idx][4] = itype;
            X->ParaPairExcitationOperator[idx] = dvalue_re + I * dvalue_im;
          }
          else{
            X->PairExcitationOperator[idx][0] = isite2;
            X->PairExcitationOperator[idx][1] = isigma2;
            X->PairExcitationOperator[idx][2] = isite1;
            X->PairExcitationOperator[idx][3] = isigma1;
            X->PairExcitationOperator[idx][4] = itype;
            X->ParaPairExcitationOperator[idx] = -(dvalue_re + I * dvalue_im);
          }
          idx++;
        }
        if (idx != X->NPairExcitationOperator) {
          fclose(fp);
          return ReadDefFileError(defname);
        }
      }
      break;

    default:
      break;
    }
    fclose(fp);

    switch(iKWidx){
    case KWCoulombIntra:
    case KWCoulombInter:
    case KWHund:
    case KWPairHop:
    case KWExchange:
    case KWIsing:
    case KWPairLift:
      if(X->iFlgGeneralSpin==TRUE){
        fprintf(stdoutMPI, "%s", cErrIncorrectFormatInter);
        return(-1);
      }
      break;
    default:
      break;
    }
  }

  ResetInteractionNum(X);
  /*=======================================================================*/
  return 0;
}

/**
 * @brief Check Site Number.
 * @param[in] *iSite a site number.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckSite(
              const int iSite,
              const int iMaxNum
              )
{
  if(iSite>=iMaxNum) return(-1);
  return 0;
}

/**
 * @brief Check Site Number for a pair -> (siteA, siteB).
 * @param[in] iSite1 a site number on a site A.
 * @param[in] iSite2 a site number on a site B.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckPairSite(
                  const int iSite1,
                  const int iSite2,
                  const int iMaxNum
                  )
{
  if(CheckSite(iSite1, iMaxNum)!=0){
    return(-1);
  }
  if(CheckSite(iSite2, iMaxNum)!=0){
    return(-1);
  }
  return 0;
}

/**
 * @brief Check Site Number for a quad -> (siteA, siteB, siteC, siteD).
 * @param[in] iSite1 a site number on site A.
 * @param[in] iSite2 a site number on site B.
 * @param[in] iSite3 a site number on site C.
 * @param[in] iSite4 a site number on site D.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckQuadSite(
                  const int iSite1,
                  const int iSite2,
                  const int iSite3,
                  const int iSite4,
                  const int iMaxNum
                  )
{
  if(CheckPairSite(iSite1, iSite2, iMaxNum)!=0){
    return(-1);
  }
  if(CheckPairSite(iSite3, iSite4, iMaxNum)!=0){
    return(-1);
  }
  return 0;
}

/**
 * @brief Check Hermite for Transfer integrals.
 * @param[in] X Define List for getting transfer integrals.
 * @retval 0 Hermite.
 * @retval -1 NonHermite.
 * @version 0.2
 * @details rearray a GeneralTransfer array to satisfy a condition of hermite conjugation between 2*i and 2*i+1 components.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckTransferHermite
(
 struct DefineList *X

 )
{
  unsigned int i,j;
  int isite1, isite2;
  int isigma1, isigma2;
  int itmpsite1, itmpsite2;
  int itmpsigma1, itmpsigma2;
  int itmperrsite1, itmperrsite2;
  int itmperrsigma1, itmperrsigma2;
  double complex dcerrTrans;
  int icheckHermiteCount=FALSE;
  int iCount=0;

  double  complex ddiff_trans;
  unsigned int itmpIdx, icntHermite, icntchemi;
  icntHermite=0;
  icntchemi=0;
  
  for(i=0; i<X->NTransfer; i++){
    isite1=X->GeneralTransfer[i][0];
    isigma1=X->GeneralTransfer[i][1];
    isite2=X->GeneralTransfer[i][2];
    isigma2=X->GeneralTransfer[i][3];
    icheckHermiteCount=FALSE;
   // fprintf(stdoutMPI, "Debug: isite1=%d, isigma1=%d, isite2=%d, isigma2=%d, reTrans=%lf, imTrans = %lf\n",
   //         isite1, isigma1, isite2, isigma2, creal(X->ParaGeneralTransfer[i]), cimag((X->ParaGeneralTransfer[i])));
    for(j=0; j<X->NTransfer; j++){
      itmpsite1=X->GeneralTransfer[j][0];
      itmpsigma1=X->GeneralTransfer[j][1];
      itmpsite2=X->GeneralTransfer[j][2];
      itmpsigma2=X->GeneralTransfer[j][3];
      if(isite1 == itmpsite2 && isite2 == itmpsite1){
        if(isigma1 == itmpsigma2 && isigma2 == itmpsigma1){
	  
          ddiff_trans = X->ParaGeneralTransfer[i]-conj(X->ParaGeneralTransfer[j]);
          if(cabs(ddiff_trans) > eps_CheckImag0 ){
            itmperrsite1=itmpsite1;
            itmperrsigma1=itmpsigma1;
            itmperrsite2=itmpsite2;
            itmperrsigma2=itmpsigma2;
            dcerrTrans=X->ParaGeneralTransfer[j];
            fprintf(stdoutMPI, cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, creal(X->ParaGeneralTransfer[i]), cimag(X->ParaGeneralTransfer[i]));
            fprintf(stdoutMPI, cErrNonHermiteTrans, itmperrsite1, itmperrsigma1, itmperrsite2, itmperrsigma2, creal(dcerrTrans), cimag(dcerrTrans));
            iCount++;
          }
          else{
            if (icheckHermiteCount == FALSE) {	      
              if(i<=j){
                if(2*icntHermite >= X->NTransfer){
                  fprintf(stderr, "Elements of Transfers are incorrect.\n");
                  return(-1);
                }
                if(isite1 !=isite2 || isigma1 !=isigma2){
                  for(itmpIdx=0; itmpIdx<4; itmpIdx++){
                    X->EDGeneralTransfer[2*icntHermite][itmpIdx]=X->GeneralTransfer[i][itmpIdx];
                    X->EDGeneralTransfer[2*icntHermite+1][itmpIdx]=X->GeneralTransfer[j][itmpIdx];
                  }
                  X->EDParaGeneralTransfer[2*icntHermite]=X->ParaGeneralTransfer[i];
                  X->EDParaGeneralTransfer[2*icntHermite+1]=X->ParaGeneralTransfer[j];
                  icntHermite++;
                }
                else{
                  X->EDChemi[icntchemi]     = X->GeneralTransfer[i][0];      
                  X->EDSpinChemi[icntchemi] = X->GeneralTransfer[i][1];      
                  X->EDParaChemi[icntchemi] = creal(X->ParaGeneralTransfer[i]);
                  icntchemi+=1;
                }
              }
              icheckHermiteCount = TRUE;
            }
          }
        }  
      
      }
    }

    //if counterpart for satisfying hermite conjugate does not exist.
    if(icheckHermiteCount == FALSE){
      fprintf(stdoutMPI, cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, creal(X->ParaGeneralTransfer[i]), cimag(X->ParaGeneralTransfer[i]));
      iCount++;
    }
  }

  if(iCount !=0){
    return -1;
  }
  X->EDNTransfer=2*icntHermite;
  X->EDNChemi=icntchemi;

  /*
    for(i=0; i<X->EDNTransfer; i++){
    for(itmpIdx=0; itmpIdx<4; itmpIdx++){
    X->GeneralTransfer[i][itmpIdx]=X->EDGeneralTransfer[i][itmpIdx];
    }
    X->ParaGeneralTransfer[i]=X->EDParaGeneralTransfer[i];
    } 
  */

  return 0;
}

/** 
 * @brief function of checking hermite conditions about interall interactions
 * 
 * @param X define list to get interall off diagonal interactions
 * 
 * @retval 0 Hermite condition is satisfied
 * @retval -1 Hermite condition is not satisfied
 * @version 0.2
 * @details rearray a InterAll_OffDiagonal array to satisfy a condition of hermite conjugation between 2*i and 2*i+1 components.
 * 
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckInterAllHermite
(
 const struct DefineList *X
 ){
  unsigned int i,j, icntincorrect, itmpret;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  int itmpsite1, itmpsite2, itmpsite3, itmpsite4;
  int itmpsigma1, itmpsigma2, itmpsigma3, itmpsigma4;
  unsigned int itmpIdx, icntHermite;
  int icheckHermiteCount=FALSE;
  double  complex ddiff_intall;
  icntincorrect=0;
  icntHermite=0;
  for(i=0; i<X->NInterAll_OffDiagonal; i++){
    itmpret=0;
    isite1=X->InterAll_OffDiagonal[i][0];
    isigma1=X->InterAll_OffDiagonal[i][1];
    isite2=X->InterAll_OffDiagonal[i][2];
    isigma2=X->InterAll_OffDiagonal[i][3];
    isite3=X->InterAll_OffDiagonal[i][4];
    isigma3=X->InterAll_OffDiagonal[i][5];
    isite4=X->InterAll_OffDiagonal[i][6];
    isigma4=X->InterAll_OffDiagonal[i][7];
    icheckHermiteCount=FALSE;

    for(j=0; j<X->NInterAll_OffDiagonal; j++) {
      itmpsite1 = X->InterAll_OffDiagonal[j][0];
      itmpsigma1 = X->InterAll_OffDiagonal[j][1];
      itmpsite2 = X->InterAll_OffDiagonal[j][2];
      itmpsigma2 = X->InterAll_OffDiagonal[j][3];
      itmpsite3 = X->InterAll_OffDiagonal[j][4];
      itmpsigma3 = X->InterAll_OffDiagonal[j][5];
      itmpsite4 = X->InterAll_OffDiagonal[j][6];
      itmpsigma4 = X->InterAll_OffDiagonal[j][7];
      
      if (isite1 == itmpsite4 && isite2 == itmpsite3 && isite3 == itmpsite2 && isite4 == itmpsite1) {
        if (isigma1 == itmpsigma4 && isigma2 == itmpsigma3 && isigma3 == itmpsigma2 && isigma4 == itmpsigma1) {
          ddiff_intall = cabs(X->ParaInterAll_OffDiagonal[i] - conj(X->ParaInterAll_OffDiagonal[j]));
	  
          if (cabs(ddiff_intall) < eps_CheckImag0) {
            itmpret=1;
            if (icheckHermiteCount == FALSE) {
              icheckHermiteCount = TRUE; //for not double counting
              if (i <= j) {
                if (2 * icntHermite >= X->NInterAll_OffDiagonal) {
                  fprintf(stdoutMPI, "Elements of InterAll are incorrect.\n");
                  return(-1);
                }
		
                for (itmpIdx = 0; itmpIdx < 8; itmpIdx++) {
                  X->InterAll[2 * icntHermite][itmpIdx] = X->InterAll_OffDiagonal[i][itmpIdx];
                  X->InterAll[2 * icntHermite + 1][itmpIdx] = X->InterAll_OffDiagonal[j][itmpIdx];
                }
		
                X->ParaInterAll[2 * icntHermite] = X->ParaInterAll_OffDiagonal[i];
                X->ParaInterAll[2 * icntHermite + 1] = X->ParaInterAll_OffDiagonal[j];
                icntHermite++;
              }
              break;
            }
          }
        }
      }
      else if (isite1 == itmpsite2 && isite2 == itmpsite1 && isite3 == itmpsite4 && isite4 == itmpsite3) {      //for spin and Kondo
          if(X->iCalcModel == Kondo ||X->iCalcModel == KondoGC || X->iCalcModel == Spin || X->iCalcModel == SpinGC) {
              if (isigma1 == itmpsigma2 && isigma2 == itmpsigma1 && isigma3 == itmpsigma4 && isigma4 == itmpsigma3) {
                  ddiff_intall = X->ParaInterAll_OffDiagonal[i] - conj(X->ParaInterAll_OffDiagonal[j]);
                  if (cabs(ddiff_intall) < eps_CheckImag0) {
                      itmpret = 1;
                      if (icheckHermiteCount == FALSE) {
                          icheckHermiteCount = TRUE; // for not double-counting
                          if (i <= j) {
                              if (2 * icntHermite >= X->NInterAll_OffDiagonal) {
                                  fprintf(stdoutMPI, "Elements of InterAll are incorrect.\n");
                                  return (-1);
                              }
                              for (itmpIdx = 0; itmpIdx < 8; itmpIdx++) {
                                  X->InterAll[2 * icntHermite][itmpIdx] = X->InterAll_OffDiagonal[i][itmpIdx];
                              }
                              for (itmpIdx = 0; itmpIdx < 4; itmpIdx++) {
                                  X->InterAll[2 * icntHermite + 1][2 * itmpIdx] = X->InterAll_OffDiagonal[i][6 -
                                                                                                             2 *
                                                                                                             itmpIdx];
                                  X->InterAll[2 * icntHermite + 1][2 * itmpIdx + 1] = X->InterAll_OffDiagonal[i][7 - 2 *
                                                                                                                     itmpIdx];

                              }
                              X->ParaInterAll[2 * icntHermite] = X->ParaInterAll_OffDiagonal[i];
                              X->ParaInterAll[2 * icntHermite + 1] = X->ParaInterAll_OffDiagonal[j];
                              icntHermite++;
                          }
                          break;
                      }
                  }
              }
          }
      }
    }
    
    //if counterpart for satisfying hermite conjugate does not exist.
    if(itmpret !=1){
      fprintf(stdoutMPI, cErrNonHermiteInterAll, isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4, creal(X->ParaInterAll_OffDiagonal[i]), cimag(X->ParaInterAll_OffDiagonal[i]));
      icntincorrect++;
    }    
  }
  
  if(icntincorrect !=0 ){
    return(-1);
  }

  /* for debug
     if((2*icntHermite) != X->NInterAll_OffDiagonal){  
     return(-1);
     }
  */
    
  for(i=0; i<X->NInterAll_OffDiagonal; i++){
    for(itmpIdx=0; itmpIdx<8; itmpIdx++){
      X->InterAll_OffDiagonal[i][itmpIdx]=X->InterAll[i][itmpIdx];
    }
    X->ParaInterAll_OffDiagonal[i]=X->ParaInterAll[i];
  }
    
  return 0;
}



/** 
 * @brief function of getting diagonal components form interall interactions
 * 
 * @param[in] X define list to get information of interall interactions
 * 
 * @retval 0  succeed to get diagonal interactions
 * @retval -1 format of interall interactions is incorrect
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
/*
int GetDiagonalInterAll
(
 struct DefineList *X
 )
{
  unsigned int i,icnt_diagonal, icnt_offdiagonal, tmp_i;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  int iret=0;
  icnt_diagonal=0;
  icnt_offdiagonal=0;
  
  setmem_IntAll_Diagonal(X);
  for(i=0; i<X->NInterAll; i++){
    isite1=X->InterAll[i][0];
    isigma1=X->InterAll[i][1];
    isite2=X->InterAll[i][2];
    isigma2=X->InterAll[i][3];
    isite3=X->InterAll[i][4];
    isigma3=X->InterAll[i][5];
    isite4=X->InterAll[i][6];
    isigma4=X->InterAll[i][7];
    
    //Get Diagonal term
    if(isite1 == isite2 && isite3 == isite4 &&
       isigma1 == isigma2  && isigma3 == isigma4)
      {
        X->InterAll_Diagonal[icnt_diagonal][0]=isite1;
        X->InterAll_Diagonal[icnt_diagonal][1]=isigma1;
        X->InterAll_Diagonal[icnt_diagonal][2]=isite3;
        X->InterAll_Diagonal[icnt_diagonal][3]=isigma3;
        X->ParaInterAll_Diagonal[icnt_diagonal] = creal(X->ParaInterAll[i]);
        icnt_diagonal++;
        continue;
      }
    else if(isite1 == isite4 && isite2 ==isite3 &&
            isigma1 == isigma4 && isigma2 ==isigma3)
      {
        X->InterAll_Diagonal[icnt_diagonal][0]=isite1;
        X->InterAll_Diagonal[icnt_diagonal][1]=isigma1;
        X->InterAll_Diagonal[icnt_diagonal][2]=isite2;
        X->InterAll_Diagonal[icnt_diagonal][3]=isigma2;
        X->ParaInterAll_Diagonal[icnt_diagonal] = -creal(X->ParaInterAll[i]);
        X->EDChemi[X->EDNChemi]     = isite1;      
        X->EDSpinChemi[X->EDNChemi] = isigma1;
        //transfer integral has minus sign for default setting
        X->EDParaChemi[X->EDNChemi] = -creal(X->ParaInterAll[i]); 
        icnt_diagonal++;
        X->EDNChemi +=1;
        continue;
      }
    else{
      //Get Off-Diagonal term
      switch(X->iCalcModel){
      case Hubbard:
      case HubbardNConserved:
      case Kondo:
      case KondoGC:
      case HubbardGC:
        if(isigma1 == isigma2 && isigma3 == isigma4){
          for(tmp_i=0; tmp_i<8; tmp_i++){
            X->InterAll_OffDiagonal[icnt_offdiagonal][tmp_i]=X->InterAll[i][tmp_i];
          }
          X->ParaInterAll_OffDiagonal[icnt_offdiagonal] = X->ParaInterAll[i];
        }
        else if(isigma1==isigma4 && isigma2 == isigma3){
          X->InterAll_OffDiagonal[icnt_offdiagonal][0]=isite1;
          X->InterAll_OffDiagonal[icnt_offdiagonal][1]=isigma1;
          X->InterAll_OffDiagonal[icnt_offdiagonal][2]=isite4;
          X->InterAll_OffDiagonal[icnt_offdiagonal][3]=isigma1;
          X->InterAll_OffDiagonal[icnt_offdiagonal][4]=isite3;
          X->InterAll_OffDiagonal[icnt_offdiagonal][5]=isigma2;
          X->InterAll_OffDiagonal[icnt_offdiagonal][6]=isite2;
          X->InterAll_OffDiagonal[icnt_offdiagonal][7]=isigma2;
          X->ParaInterAll_OffDiagonal[icnt_offdiagonal] = -X->ParaInterAll[i];
        }
        else{
          // Sz symmetry is assumed
          if(X->iCalcModel==Hubbard || X->iCalcModel==Kondo){
            fprintf(stdoutMPI, cErrNonConservedInterAll,
                    isite1,
                    isigma1,
                    isite2,
                    isigma2,
                    isite3,
                    isigma3,
                    isite4,
                    isigma4,
                    creal(X->ParaInterAll[i]),
                    cimag(X->ParaInterAll[i])
                    );
            iret=-1;
          }
          else{
            for(tmp_i=0; tmp_i<8; tmp_i++){
              X->InterAll_OffDiagonal[icnt_offdiagonal][tmp_i]=X->InterAll[i][tmp_i];
            }
            X->ParaInterAll_OffDiagonal[icnt_offdiagonal] = X->ParaInterAll[i];
          }
        }
        break;
      case Spin:
      case SpinGC:
        if(isite1 == isite2 && isite3 == isite4){
          for(tmp_i=0; tmp_i<8; tmp_i++){
            X->InterAll_OffDiagonal[icnt_offdiagonal][tmp_i]=X->InterAll[i][tmp_i];
          }
          X->ParaInterAll_OffDiagonal[icnt_offdiagonal] = X->ParaInterAll[i];
        }      
        break;
      default:
        return(-1);
      }
      if(iret != -1){
        icnt_offdiagonal++;
      }
    }

    if(iret !=0){
      return(-1);
    }
  }
  
  return 0;
}
*/

/**
 * @brief function of getting diagonal components form Time-dependent interall interactions
 *
 * @param[in] X define list to get information of Time-dependent interall interactions
 *
 * @retval 0  succeed to get diagonal interactions
 * @retval -1 format of interall interactions is incorrect
 * @version 2.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int GetDiagonalInterAll
        (
                int **InterAll,
                complex double *ParaInterAll,
                const int NInterAll,
                int **InterAllDiagonal,
                double *ParaInterAllDiagonal,
                int **InterAllOffDiagonal,
                complex double *ParaInterAllOffDiagonal,
                int *Chemi,
                int *SpinChemi,
                double *ParaChemi,
                int *NChemi,
                const int iCalcModel
        )
{
  unsigned int i,icnt_diagonal, icnt_offdiagonal, tmp_i;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  int iret=0;
  icnt_diagonal=0;
  icnt_offdiagonal=0;

  setmem_IntAll_Diagonal(
          InterAllOffDiagonal, ParaInterAllOffDiagonal,
          InterAllDiagonal, ParaInterAllDiagonal, NInterAll);

  for(i=0; i<NInterAll; i++){
    isite1=InterAll[i][0];
    isigma1=InterAll[i][1];
    isite2=InterAll[i][2];
    isigma2=InterAll[i][3];
    isite3=InterAll[i][4];
    isigma3=InterAll[i][5];
    isite4=InterAll[i][6];
    isigma4=InterAll[i][7];

    //Get Diagonal term
    if(isite1 == isite2 && isite3 == isite4 &&
       isigma1 == isigma2  && isigma3 == isigma4)
    {
      InterAllDiagonal[icnt_diagonal][0]=isite1;
      InterAllDiagonal[icnt_diagonal][1]=isigma1;
      InterAllDiagonal[icnt_diagonal][2]=isite3;
      InterAllDiagonal[icnt_diagonal][3]=isigma3;
      ParaInterAllDiagonal[icnt_diagonal] = creal(ParaInterAll[i]);
      icnt_diagonal++;
      continue;
    }
    else if(isite1 == isite4 && isite2 ==isite3 &&
            isigma1 == isigma4 && isigma2 ==isigma3)
    {
      InterAllDiagonal[icnt_diagonal][0]=isite1;
      InterAllDiagonal[icnt_diagonal][1]=isigma1;
      InterAllDiagonal[icnt_diagonal][2]=isite2;
      InterAllDiagonal[icnt_diagonal][3]=isigma2;
      ParaInterAllDiagonal[icnt_diagonal] = -creal(ParaInterAll[i]);
      Chemi[*NChemi]     = isite1;
      SpinChemi[*NChemi] = isigma1;
      //transfer integral has minus sign for default setting
      ParaChemi[*NChemi] = -creal(ParaInterAll[i]);
      icnt_diagonal++;
      *NChemi +=1;
      continue;
    }
    else{
      //Get Off-Diagonal term
      switch(iCalcModel){
        case Hubbard:
        case HubbardNConserved:
        case Kondo:
        case KondoGC:
        case HubbardGC:
          if(isigma1 == isigma2 && isigma3 == isigma4){
            for(tmp_i=0; tmp_i<8; tmp_i++){
              InterAllOffDiagonal[icnt_offdiagonal][tmp_i]=InterAll[i][tmp_i];
            }
            ParaInterAllOffDiagonal[icnt_offdiagonal] = ParaInterAll[i];
          }
          else if(isigma1==isigma4 && isigma2 == isigma3){
            InterAllOffDiagonal[icnt_offdiagonal][0]=isite1;
            InterAllOffDiagonal[icnt_offdiagonal][1]=isigma1;
            InterAllOffDiagonal[icnt_offdiagonal][2]=isite4;
            InterAllOffDiagonal[icnt_offdiagonal][3]=isigma1;
            InterAllOffDiagonal[icnt_offdiagonal][4]=isite3;
            InterAllOffDiagonal[icnt_offdiagonal][5]=isigma2;
            InterAllOffDiagonal[icnt_offdiagonal][6]=isite2;
            InterAllOffDiagonal[icnt_offdiagonal][7]=isigma2;
            ParaInterAllOffDiagonal[icnt_offdiagonal] = -ParaInterAll[i];
          }
          else{
            // Sz symmetry is assumed
            if(iCalcModel==Hubbard || iCalcModel==Kondo){
              fprintf(stdoutMPI, cErrNonConservedInterAll,
                      isite1,
                      isigma1,
                      isite2,
                      isigma2,
                      isite3,
                      isigma3,
                      isite4,
                      isigma4,
                      creal(ParaInterAll[i]),
                      cimag(ParaInterAll[i])
              );
              iret=-1;
            }
            else{
              for(tmp_i=0; tmp_i<8; tmp_i++){
                InterAllOffDiagonal[icnt_offdiagonal][tmp_i]=InterAll[i][tmp_i];
              }
              ParaInterAllOffDiagonal[icnt_offdiagonal] = ParaInterAll[i];
            }
          }
          break;
        case Spin:
        case SpinGC:
          if(isite1 == isite2 && isite3 == isite4){
            for(tmp_i=0; tmp_i<8; tmp_i++){
              InterAllOffDiagonal[icnt_offdiagonal][tmp_i]=InterAll[i][tmp_i];
            }
              ParaInterAllOffDiagonal[icnt_offdiagonal] =ParaInterAll[i];
          }
          break;
        default:
          return(-1);
      }
      if(iret != -1){
        icnt_offdiagonal++;
      }
    }

    if(iret !=0){
      return(-1);
    }
  }

  return 0;
}


/** 
 * @brief function of judging a type of define files.
 * 
 * @param[in] argc argument count
 * @param[in] argv argument vector 
 * @param[out] mode a number to show a type of a define file
 * 
 * @retval 0 format is correct
 * @retval -1 format is incorrect
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int JudgeDefType
(
 const int argc,
 char *argv[],
 int *mode
 )
{
   int ver_maj =
#include "version_major.h"
;
   int ver_min =
#include "version_miner.h"
;
   int ver_pat =
#include "version_patch.h"
;

  if(argc == 3 && 
     (CheckWords(argv[1], "-e") == 0 ||
      CheckWords(argv[1], "--Expert") == 0)){
    *mode=EXPERT_MODE;
  }
  else if (argc == 3 && 
           (CheckWords(argv[1], "-s") ==0 ||
            CheckWords(argv[1], "--Standard") == 0 )){
    *mode=STANDARD_MODE;
  }
  else if (argc == 3 && 
           (CheckWords(argv[1], "-sdry") == 0 ||
            CheckWords(argv[1], "-s-dry") == 0)
           ){
    *mode = STANDARD_DRY_MODE;
  }
  else if (argc >= 2 &&
           (CheckWords(argv[1], "-v") == 0
            || CheckWords(argv[1], "--version") == 0)
           ) {
    fprintf(stdoutMPI, "\nHPhi version %d.%d.%d \n\n", ver_maj, ver_min, ver_pat);
    exit(-1);
  }
  else{
    /*fprintf(stdoutMPI, cErrArgv, argv[1]);*/
    fprintf(stdoutMPI, "\n[Usage] \n");
    fprintf(stdoutMPI, "* Expert mode \n");
    fprintf(stdoutMPI, "   $ HPhi -e {namelist_file} \n");
    fprintf(stdoutMPI, "* Standard mode \n");
    fprintf(stdoutMPI, "   $ HPhi -s {input_file} \n");
    fprintf(stdoutMPI, "* Standard DRY mode \n");
    fprintf(stdoutMPI, "   $ HPhi -sdry {input_file} \n");
    fprintf(stdoutMPI, "   In this mode, Hphi stops after it generats expert input files. \n");
    fprintf(stdoutMPI, "* Print the version \n");
    fprintf(stdoutMPI, "   $ HPhi -v \n\n");
    exit(-1);
  }

  return 0;
}

/** 
 * @brief function of checking format of spin interactions
 * 
 * @param[in] site1 a site number on site1.
 * @param[in] site2 a site number on site2.
 * @param[in] site3 a site number on site3.
 * @param[in] site4 a site number on site4.
 * 
 * @retval 0 format is correct
 * @retval -1 format is incorrect
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckFormatForSpinInt
(
 const int site1,
 const int site2,
 const int site3,
 const int site4
 ){
  if(site1==site2 && site3==site4){
    return 0;
  }

  fprintf(stdoutMPI, cWarningIncorrectFormatForSpin, site1, site2, site3, site4);
  return(-1);

}

/** 
 * 
 * @brief function of checking format of Kondo interactions
 * 
 * @param[in] X define list to get information of interall interactions
 * 
 * @retval 0 format is correct
 * @retval -1 format is incorrect
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckFormatForKondoInt
        (
                const int isite1, const int isite2,
                const int isite3, const int isite4,
                int* iLocInfo
        )
{
  if (iLocInfo[isite1] != ITINERANT || iLocInfo[isite2] != ITINERANT) {
    if (isite1 != isite2) {
      fprintf(stdoutMPI, cErrIncorrectFormatForKondoInt, isite1, isite2, isite3, isite4);
      return -1;
    }
  }
  if (iLocInfo[isite3] != ITINERANT || iLocInfo[isite4] != ITINERANT) {
    if (isite3 != isite4) {
      fprintf(stdoutMPI, cErrIncorrectFormatForKondoInt, isite1, isite2, isite3, isite4);
      return -1;
    }
  }
  return 0;
}

/*
int CheckFormatForKondoInt
(
 struct DefineList *X
 ){
  unsigned int i,iboolLoc;
  int isite1, isite2, isite3, isite4;

  iboolLoc=0;
  for(i=0; i<X->NInterAll; i++){
    isite1=X->InterAll[i][0];
    isite2=X->InterAll[i][2];
    isite3=X->InterAll[i][4];
    isite4=X->InterAll[i][6];

    if(X->LocSpn[isite1]!=ITINERANT || X->LocSpn[isite2]!=ITINERANT){
      if(isite1 != isite2){
        iboolLoc=1;
        fprintf(stdoutMPI, cErrIncorrectFormatForKondoInt, isite1, isite2, isite3, isite4);
        continue;
      }
    }
    if(X->LocSpn[isite3]!=ITINERANT || X->LocSpn[isite4]!=ITINERANT){
      if(isite3 != isite4){
        iboolLoc=1;
        fprintf(stdoutMPI, cErrIncorrectFormatForKondoInt, isite1, isite2, isite3, isite4);
        continue;
      }
    }
  }
  if(iboolLoc==1){
    return(-1);
  }  
  return 0;
}
*/

/** 
 * @brief function to set convergence factors
 * 
 * @param[in] X Define list to get Lanczos eps.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void SetConvergenceFactor
(
 struct DefineList *X
 )
{
  //In future, convergence facator can be set by a def file.
  int neps = -8;
  int nepsCG =-8;
  int nbisec =-14;
  int nEnergy = -12;
  int nShiftBeta=8;
  int nepsvec12=-14;
  eps=pow(10.0, neps);
  eps_CG=pow(10.0, nepsCG);
  eps_Lanczos     = pow(10,-X->LanczosEps);
  eps_Bisec = pow(10.0, nbisec);
  eps_Energy = pow(10.0, nEnergy);
  dShiftBeta = pow(10.0, nShiftBeta);
  eps_vec12 = pow(10.0, nepsvec12);
}

/** 
 * @brief function of checking indexies of localized spin
 * 
 * @param[in/out] X Define list to get and put information of localized spin
 * 
 * @return TURE Indecies of localizes spin is correct
 * @return FALSE Indecies of localizes spin is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckLocSpin
(
 struct DefineList *X
 )
{

  unsigned int i=0;
  switch(X->iCalcModel){
  case Hubbard:
  case HubbardNConserved:
  case HubbardGC:
  case SpinlessFermion:
  case SpinlessFermionGC:
    for(i=0; i<X->Nsite; i++){
      if(X->LocSpn[i]!=ITINERANT){
        return FALSE;
      }
    }
    break;

  case Kondo:
  case KondoGC:
    for(i=0; i<X->Nsite; i++){
      if(X->LocSpn[i]>LOCSPIN){
        X->iFlgGeneralSpin=TRUE;
      }
      else if(X->LocSpn[i]<ITINERANT){
        return FALSE;
      }
    }
    break;

  case Spin:
  case SpinGC:
    for(i=0; i<X->Nsite; i++){
      if(X->LocSpn[i]>LOCSPIN){
        X->iFlgGeneralSpin=TRUE;
      }
      else if(X->LocSpn[i]<LOCSPIN){
        return FALSE;
      }
    }
    break;
  
  default:
    return FALSE;
    //break;
  }

  if(CheckTotal2Sz(X) != TRUE){
    return FALSE;
  }
  return TRUE;
}  

/** 
 * 
 * @brief function of resetting number of interactions
 * 
 * @param[out] X Define list to add number of ising coulomnb interactions
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
void ResetInteractionNum
(
 struct DefineList *X
 )
{
  X->NHundCoupling += X->NIsingCoupling;
  X->NCoulombInter += X->NIsingCoupling;
}

/** 
 * @brief function of initializeing interactions
 * 
 * @param[out] X Define list to initialize number of interactions
 * @version 0.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
void InitializeInteractionNum
(
 struct DefineList *X
 )
{
  X->NTransfer=0;
  X->NCoulombIntra=0;
  X->NCoulombInter=0;
  X->NIsingCoupling=0;
  X->NPairLiftCoupling=0;
  X->NInterAll=0;
  X->NCisAjt=0;
  X->NCisAjtCkuAlvDC=0;
  X->NSingleExcitationOperator=0;
  X->NPairExcitationOperator=0;
  //[s] Time Evolution
  X->NTETimeSteps=0;
  X->NLaser=0;
  X->NTEInterAll=0;
  X->NTETransfer=0;
  //[e] Time Evolution

}

/** 
 * @brief function of checking spin index for all interactions
 * 
 * @param[in] X Define list to get informations of all interactions
 * @retval TRUE spin index is correct
 * @retval FALSE spin index is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckGeneralSpinIndexForInterAll
(
        const int isite1, const int isigma1,
        const int isite2, const int isigma2,
        const int isite3, const int isigma3,
        const int isite4, const int isigma4,
        int* iLocInfo
 )
{
   if( isigma1 > iLocInfo[isite1] || isigma2 >iLocInfo[isite2]
         ||isigma3 > iLocInfo[isite3] || isigma4 >iLocInfo[isite4]){
        fprintf(stdoutMPI, "%s", cErrIncorrectSpinIndexForInter);
        return FALSE;
    }
  return TRUE;
}

/** 
 * @brief function of checking spin index for transfers
 * 
 * @param[in] X Define list to get informations of transfers
 * @retval TRUE spin index is correct
 * @retval FALSE spin index is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckSpinIndexForTrans
(
 struct DefineList *X
 )
{
  unsigned int i=0;
  int isite1, isite2;
  int isigma1, isigma2;
  if(X->iFlgGeneralSpin==TRUE){
    for(i=0; i<X->NTransfer; i++){
      isite1 =X->GeneralTransfer[i][0];
      isigma1=X->GeneralTransfer[i][1];
      isite2 =X->GeneralTransfer[i][2];
      isigma2=X->GeneralTransfer[i][3];
      if(isigma1 > X->LocSpn[isite1] || isigma2 >X->LocSpn[isite2]){
        fprintf(stdoutMPI, "%s", cErrIncorrectSpinIndexForTrans);
        return FALSE;
      }
    }
  }
  return TRUE;
}

/** 
 * @brief function of checking an input data of total2Sz
 * 
 * @param[in] X Define list to get informations of transfers
 * @retval TRUE spin index is correct
 * @retval FALSE spin index is incorrect
 * @version 0.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Takahiro Misawa (The University of Tokyo)
 */
int CheckTotal2Sz
(
 struct DefineList *X
 )
{
  if(X->iFlgSzConserved==TRUE && X->iFlgGeneralSpin==FALSE){
    int tmp_Nup=X->NLocSpn+X->NCond+X->Total2Sz;
    int tmp_Ndown=X->NLocSpn+X->NCond-X->Total2Sz;
    if(tmp_Nup%2 != 0 && tmp_Ndown%2 !=0){
      printf("Nup=%d, Ndown=%d\n",X->Nup,X->Ndown);
      fprintf(stdoutMPI, "2Sz is incorrect.\n");
      return FALSE;
    }
  }
  return TRUE;
}

/** 
 * 
 * @brief function of checking whether ctmp is same as cKeyWord or not
 * 
 * @param[in] ctmp 
 * @param[in] cKeyWord 
 * @return 0 ctmp is same as cKeyWord
 * 
 * @version 1.1.0
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckWords(
               const char* ctmp,
               const char* cKeyWord
               )
{
  unsigned int i=0;

  char ctmp_small[256]={0};
  char cKW_small[256]={0};
  unsigned int n;
  n=strlen(cKeyWord);
  strncpy(cKW_small, cKeyWord, n);

  for(i=0; i<n; i++){
    cKW_small[i]=tolower(cKW_small[i]);
  }
  n=strlen(ctmp);
  strncpy(ctmp_small, ctmp, n);
  for(i=0; i<n; i++){
    ctmp_small[i]=tolower(ctmp_small[i]);
  }
  if(n<strlen(cKW_small)) n=strlen(cKW_small);
  return(strncmp(ctmp_small, cKW_small, n));
}

int GetFileNameByKW(
        int iKWidx,
        char **FileName
){
  if(cFileNameListFile == NULL){
    return -1;
  }
  *FileName=cFileNameListFile[iKWidx];
  return 0;
}


/**
 * @brief Check InterAll condition.
 * @param[in] iCalcModel Target Model defined in CalcMod file (ex. Spin, SpinGC etc.).
 * @param[in] Nsite  A total number of site.
 * @param[in] isite1 a site number on the site A.
 * @param[in] isigma1 a spin index on the site A.
 * @param[in] isite2 a site number on the site B.
 * @param[in] isigma2 a spin index on the site B.
 * @param[in] isite3 a site number on the site C.
 * @param[in] isigma3 a spin index on the site C.
 * @param[in] isite4 a site number on the site D.
 * @param[in] isigma4 a spin index on the site D.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 2.2
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckInterAllCondition(
        int iCalcModel,
        int Nsite,
        int iFlgGeneralSpin,
        int *iLocInfo,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4
){
  if(CheckQuadSite(isite1, isite2, isite3, isite4, Nsite) !=0){
    fprintf(stderr, "%s", "Error: Site index of InterAll is incorrect.\n");
    return(-1);
  }

  if(iCalcModel == Spin || iCalcModel ==SpinGC){
    if(!CheckFormatForSpinInt(isite1, isite2, isite3, isite4)==0){
      fprintf(stderr, "%s", "Error: Spin index of InterAll is incorrect.\n");
      return(-1);
    }
  }
  else if(iCalcModel == SpinlessFermion || iCalcModel==SpinlessFermionGC){
    if(isigma1 !=0 || isigma2 != 0 || isigma3 != 0 || isigma4 !=0){
      fprintf(stderr, "%s", "Error: Spin index of InterAll is incorrect.\n");
      return -1;
    }
  }
  else if(iCalcModel == Kondo){
    if(CheckFormatForKondoInt(isite1, isite2, isite3, isite4, iLocInfo)!=0){
      return -1;
    }
  }

  if(iFlgGeneralSpin ==TRUE) {
    if(CheckGeneralSpinIndexForInterAll(isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4, iLocInfo)!=0){
      return -1;
    }
  }
  return 0;
}

int InputInterAllInfo(
        int *icnt_interall,
        int **iInterAllInfo,
        double complex *cInterAllValue,
        int isite1, int isigma1,
        int isite2, int isigma2,
        int isite3, int isigma3,
        int isite4, int isigma4,
        double dvalue_re, double dvalue_im
) {
  int i = 0;
  int iflg_interall = 0;
  //Collect and sum same components of InterAll interactions
  for (i = 0; i < *icnt_interall; i++) {
    if (isite1 == iInterAllInfo[i][0] && isite2 == iInterAllInfo[i][2] &&
        isite3 == iInterAllInfo[i][4] && isite4 == iInterAllInfo[i][6] &&
        isigma1 == iInterAllInfo[i][1] && isigma2 == iInterAllInfo[i][3] &&
        isigma3 == iInterAllInfo[i][5] && isigma4 == iInterAllInfo[i][7]) {
      cInterAllValue[i] += dvalue_re + dvalue_im * I;
      iflg_interall = 1;
      return 0;
    }
  }

  //Input all InterAll interactions
  if (iflg_interall == 0) {
    iInterAllInfo[*icnt_interall][0] = isite1;
    iInterAllInfo[*icnt_interall][1] = isigma1;
    iInterAllInfo[*icnt_interall][2] = isite2;
    iInterAllInfo[*icnt_interall][3] = isigma2;
    iInterAllInfo[*icnt_interall][4] = isite3;
    iInterAllInfo[*icnt_interall][5] = isigma3;
    iInterAllInfo[*icnt_interall][6] = isite4;
    iInterAllInfo[*icnt_interall][7] = isigma4;
    cInterAllValue[*icnt_interall] = dvalue_re + I * dvalue_im;

    //Check Diagonal part or not
    if (isite1 == isite2 && isite3 == isite4 &&
        isigma1 == isigma2 && isigma3 == isigma4) { //normal diagonal part
      return 1;
    } else if (isite1 == isite4 && isite2 == isite3 &&
               isigma1 == isigma4 && isigma2 == isigma3) { //hund term
      return 1;
    }
    *icnt_interall++;
  }
  return 0;
}



/**
 * @brief Check Hermite for TETransfer integrals.
 * @param[in] X Define List for getting transfer integrals.
 * @retval 0 Hermite.
 * @retval -1 NonHermite.
 **/
int CheckTETransferHermite
        (
                struct DefineList *X,
                const int NTETransfer,
                const int idx
        )
{
  unsigned int i,j;
  int isite1, isite2;
  int isigma1, isigma2;
  int itmpsite1, itmpsite2;
  int itmpsigma1, itmpsigma2;
  int itmperrsite1, itmperrsite2;
  int itmperrsigma1, itmperrsigma2;
  double complex dcerrTrans;
  int icheckHermiteCount;
  int iCount=0;

  double  complex ddiff_trans;
  unsigned int itmpIdx, icntHermite, icntchemi;
  icntHermite=0;
  icntchemi=0;

  int** tmp_TETransfer = (int**)malloc((NTETransfer)*sizeof(int*));
  for(i =0; i<NTETransfer; i++ ){
    tmp_TETransfer[i] = (int*)malloc((4*sizeof(int)));
  }
  double complex*tmp_paraTETransfer = (double complex*)malloc((NTETransfer)*sizeof(double complex));

  //copy
  for(i=0; i<NTETransfer; i++){
    for(j=0; j<4; j++){
      tmp_TETransfer[i][j]=X->TETransfer[idx][i][j];
      X->TETransfer[idx][i][j]=0;
    }
    tmp_paraTETransfer[i] = X->ParaTETransfer[idx][i];
    X->ParaTETransfer[idx][i]=0.0;
  }

  for(i=0; i<NTETransfer; i++){
    isite1=tmp_TETransfer[i][0];
    isigma1=tmp_TETransfer[i][1];
    isite2=tmp_TETransfer[i][2];
    isigma2=tmp_TETransfer[i][3];
    icheckHermiteCount=FALSE;
    for(j=0; j<NTETransfer; j++){
      itmpsite1=tmp_TETransfer[j][0];
      itmpsigma1=tmp_TETransfer[j][1];
      itmpsite2=tmp_TETransfer[j][2];
      itmpsigma2=tmp_TETransfer[j][3];
      if(isite1 == itmpsite2 && isite2 == itmpsite1){
        if(isigma1 == itmpsigma2 && isigma2 == itmpsigma1){

          ddiff_trans = tmp_paraTETransfer[i]-conj(tmp_paraTETransfer[j]);
          if(cabs(ddiff_trans) > eps_CheckImag0 ){
            itmperrsite1=itmpsite1;
            itmperrsigma1=itmpsigma1;
            itmperrsite2=itmpsite2;
            itmperrsigma2=itmpsigma2;
            dcerrTrans=tmp_paraTETransfer[j];
            fprintf(stdoutMPI, cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, creal(tmp_paraTETransfer[i]), cimag(tmp_paraTETransfer[i]));
            fprintf(stdoutMPI, cErrNonHermiteTrans, itmperrsite1, itmperrsigma1, itmperrsite2, itmperrsigma2, creal(dcerrTrans), cimag(dcerrTrans));
            iCount++;
          }
          else{
            if (icheckHermiteCount == FALSE) {
              if(i<=j){
                if(2*icntHermite >= NTETransfer){
                  fprintf(stderr, "Elements of Transfers are incorrect.\n");
                  return(-1);
                }
                if(isite1 !=isite2 || isigma1 !=isigma2){
                  for(itmpIdx=0; itmpIdx<4; itmpIdx++){
                    X->TETransfer[idx][2*icntHermite][itmpIdx]=tmp_TETransfer[i][itmpIdx];
                    X->TETransfer[idx][2*icntHermite+1][itmpIdx]=tmp_TETransfer[j][itmpIdx];
                  }
                  X->ParaTETransfer[idx][2*icntHermite]=tmp_paraTETransfer[i];
                  X->ParaTETransfer[idx][2*icntHermite+1]=tmp_paraTETransfer[j];
                  icntHermite++;
                }
                else{
                  X->TETransferDiagonal[idx][icntchemi][0] = tmp_TETransfer[i][0];
                  X->TETransferDiagonal[idx][icntchemi][1] = tmp_TETransfer[i][1];
                  X->ParaTETransferDiagonal[idx][icntchemi] = creal(tmp_paraTETransfer[i]);
                  icntchemi+=1;
                }
              }
              icheckHermiteCount = TRUE;
            }
          }
        }
      }
    }

    //if counterpart for satisfying hermite conjugate does not exist.
    if(icheckHermiteCount == FALSE){
      fprintf(stdoutMPI, cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, creal(tmp_paraTETransfer[i]), cimag(tmp_paraTETransfer[i]));
      iCount++;
      //fprintf(stdoutMPI, cErrNonHermiteTrans, itmperrsite1, itmperrsigma1, itmperrsite2, itmperrsigma2, creal(dcerrTrans), cimag(dcerrTrans));
      //return(-1);
    }
  }

  if(iCount !=0){
    return -1;
  }

  X->NTETransfer[idx]=2*icntHermite;
  X->NTETransferDiagonal[idx]=icntchemi;

  free(tmp_TETransfer);
  free(tmp_paraTETransfer);
  return 0;
}

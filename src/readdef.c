/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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

#include "Common.h"
#include "readdef.h"

/**
 * @brief Error Function of reading def files.
 * @param[in] *defname name of def file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileError(
		     const	char *defname
		     ){
  printf(cErrReadDefFile, defname);
  return -1;
}

/**
 * @brief Function of Validating value.
 * @param[in] icheckValue value to validate.
 * @param[in] ilowestValue lowest value which icheckValue can be set.
 * @param[in] iHighestValue heighest value which icheckValue can be set.
 * @retval 0 value is correct.
 * @retval -1 value is incorrect.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ValidateValue(
		  const int icheckValue, 
		  const int ilowestValue, 
		  const int iHighestValue
		  ){

  if(icheckValue < ilowestValue || icheckValue > iHighestValue){
    return -1;
  }
  return 0;
}

/**
 * @brief Function of Checking keyword in NameList file.
 * @param[in] cKW keyword candidate
 * @param[in] Reffercnce of keyword List
 * @param[in] number of keyword
 * @param[out] iKWidx index of keyword
 * @retval 0 keyword is correct.
 * @retval -1 keyword is incorrect.
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
    else if(strcmp(cKW, cKWList[itmpKWidx])==0){
      iret=0;
      *iKWidx=itmpKWidx;
    }
  }
  return iret;
}

/**
 * @brief Function of Getting keyword and it's variable from characters.
 * @param[in] ctmpLine characters including keyword and it's variable 
 * @param[out] keyword
 * @param[out] variable
 * @retval 0 keyword and it's variable are obtained.
 * @retval 1 ctmpLine is a comment line.
 * @retval -1 format of ctmpLine is incorrect.
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
    ctmpRead = strtok( ctmpLine, csplit);
    if(strncmp(ctmpRead, "=", 1)==0 || strncmp(ctmpRead, "#", 1)==0){
      return 1;
    }
    strcpy(ctmp, ctmpRead);
    
    ctmpRead = strtok( NULL, csplit );
    *itmp = strtol(ctmpRead, &cerror, 0);
    //if ctmpRead is not integer type
    if(*cerror != '\0'){
      printf(cErrDefFileFormat, cerror);
      return -1;
    }

    ctmpRead = strtok( NULL, csplit );
    if(ctmpRead != NULL){
      printf(cErrDefFileFormat, ctmpRead);
      return -1;
    }
    
    return 0;
}

/**
 * @brief Function of Reading calcmod file.
 * @param[in]  *defname file name to read.
 * @param[out] *X Define List for getting flags of calc-mode.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
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
  /*=======================================================================*/
  fp = fopen(defname, "r");
  if(fp==NULL) return ReadDefFileError(defname);
  /* read Parameters from calcmod.def*/
  while( fgets(ctmpLine, D_CharTmpReadDef+D_CharKWDMAX, fp)!=NULL ){
    if( (iret=GetKWWithIdx(ctmpLine, ctmp, &itmp)) !=0){
      if(iret==1) continue;
      return -1;
    }   
    if(strcmp(ctmp, "CalcType")==0){
      X->iCalcType=itmp;
    }
    else if(strcmp(ctmp, "FlgFiniteTemperature")==0){
      X->iFlgFiniteTemperature = itmp;
    }
    else if(strcmp(ctmp, "CalcModel")==0){
      X->iCalcModel=itmp;
    }
    else if(strcmp(ctmp, "OutputMode")==0){
      X->iOutputMode=itmp;
    }
    else{
      printf(cErrDefFileParam, defname, ctmp);
      return -1;
    }
  }
  fclose(fp);
  
  /* Check values*/
  if(ValidateValue(X->iCalcModel, 0, NUM_CALCMODEL-1)){
    printf(cErrCalcType, defname);
    return (-1);
  }
  if(ValidateValue(X->iCalcType, 0, NUM_CALCTYPE-1)){
    printf(cErrCalcType, defname);
    return (-1);
  }
  if(ValidateValue(X->iOutputMode, 0, NUM_OUTPUTMODE-1)){
    printf(cErrOutputMode, defname);
    return (-1);
  }
  
  /* In the case of Full Diagonalization method(iCalcType=2)*/
  if(X->iCalcType==2 && ValidateValue(X->iFlgFiniteTemperature, 0, 1)){
    printf(cErrFiniteTemp, defname);
    return (-1);
  }
  
  return 0;
}

/**
 * @brief Function of Fitting FileName
 * @param[in]  *cFileListNameFile file for getting names of input files.
 * @param[out] *cFileNameList arrays for getting names of input files.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
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
  char ctmpKW[D_CharTmpReadDef];
  int i;
  for(i=0; i< D_iKWNumDef; i++){
    strcpy(cFileNameList[i],"");
  }

  fplist = fopen(cFileListNameFile, "r");
  if(fplist==NULL) return ReadDefFileError(cFileListNameFile);

  while( fscanf(fplist,"%s %s\n", ctmpKW, ctmpFileName) !=EOF){
    if(strncmp(ctmpKW, "#", 1)==0){
      continue;
    }
    else if(strcmp(ctmpKW, "")*strcmp(ctmpFileName, "")==0){
      printf(cErrKW_InCorPair, cFileListNameFile);
      fclose(fplist);
      return -1;
    }
    /*!< Check KW */
    if( CheckKW(ctmpKW, cKWListOfFileNameList, D_iKWNumDef, &itmpKWidx)!=0 ){
      printf(cErrKW, ctmpKW, cFileListNameFile);
      printf("%s", cErrKW_ShowList);
      for(i=0; i<D_iKWNumDef;i++){
	printf("%s \n", cKWListOfFileNameList[i]);
      }
      fclose(fplist);
      return -1;
    }
    /*!< Check cFileNameList to prevent from double registering the file name */    
    if(strcmp(cFileNameList[itmpKWidx], "") !=0){
      printf(cErrKW_Same, cFileListNameFile);
      fclose(fplist);
      return -1;
    }
    /*!< Copy FileName */
    strcpy(cFileNameList[itmpKWidx], ctmpFileName);
  }
  fclose(fplist);  
  return 0;
}

/** 
 * @brief  Function of reading informations from def files.
 * @param[in] *xNameListFile List of Input File names.
 * @param[out] *X Define List for getting flags of calc-mode.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileNInt(
		    char *xNameListFile, 
		    struct DefineList *X
		    )
{
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_CharTmpReadDef];
  int itmp;
  
  printf("Start: Read File '%s'.\n", xNameListFile); 
  if(GetFileName(xNameListFile, cFileNameListFile)!=0){
    return -1;
  }
  printf("End: Read File '%s'.\n", xNameListFile); 

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
	printf(cErrMakeDef, cKWListOfFileNameList[iKWidx]);
	return -1;
	break;
      default:
	break;
      }
    }
   } 
  
  
  for(iKWidx=0; iKWidx< D_iKWNumDef; iKWidx++){ 
    strcpy(defname, cFileNameListFile[iKWidx]);

  if(strcmp(defname,"")==0) continue;
  
    printf("Read File '%s' for %s.\n", defname, cKWListOfFileNameList[iKWidx]);
    fp = fopen(defname, "r");
    if(fp==NULL) return ReadDefFileError(defname);
    switch(iKWidx){
    case KWCalcMod:
      /* Read calcmod.def---------------------------------------*/
      if(ReadcalcmodFile(defname, X)!=0){
	fclose(fp);
	 return ReadDefFileError(defname);
      }
      break;
    case KWModPara:
      /* Read modpara.def---------------------------------------*/
      //TODO: add error procedure here when parameters are not enough.
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &itmp); //2
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //3
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //4
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //5
      fscanf(fp,"%s %s\n", ctmp, X->CDataFileHead); //6
      fscanf(fp,"%s %s\n", ctmp, X->CParaFileHead); //7
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);   //8
      fscanf(fp,"%s %d\n", ctmp, &(X->Nsite));      //9
      fscanf(fp,"%s %d\n", ctmp, &(X->Nup));         //10
      fscanf(fp,"%s %d\n", ctmp, &(X->Ndown));       //11	
      if(X->iCalcModel == Spin){
	X->Ne=X->Nup;
	X->Ndown=X->Nsite-X->Nup;
      }
      fscanf(fp,"%s %d\n", ctmp, &(X->Lanczos_max)); //12
      fscanf(fp,"%s %ld\n", ctmp, &(X->initial_iv));//13
      fscanf(fp,"%s %d\n", ctmp, &(X->nvec));      //14
      fscanf(fp,"%s %d\n", ctmp, &(X->k_exct));    //15
      fscanf(fp,"%s %d\n", ctmp, &(X->LanczosEps)); //16
      fscanf(fp,"%s %d\n", ctmp, &(X->LanczosTarget)); //17
      fscanf(fp, "%s %lf\n", ctmp, &(LargeValue)); //18
      fscanf(fp, "%s %d\n", ctmp, &(NumAve)); //19
      fscanf(fp, "%s %d\n", ctmp, &(ExpecInterval)); //20
      break;
      
    case KWLocSpin:
      // Read locspn.def
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NLocSpn));
      break;
    case KWTrans: 
      // Read transfer.def
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NTransfer));
      break;
    case KWCoulombIntra:
      /* Read coulombintra.def----------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NCoulombIntra));
      break;
    case KWCoulombInter:
      /* Read coulombinter.def----------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NCoulombInter));
      break;
    case KWHund:
      /* Read hund.def------------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NHundCoupling));
      break;
    case KWPairHop:
      /* Read pairhop.def---------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NPairHopping));
      break;
    case KWExchange:
      /* Read exchange.def--------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NExchangeCoupling));
      break;
    case KWIsing:
      /* Read exchange.def--------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NIsingCoupling));
      break;
    case KWPairLift:
      /* Read exchange.def--------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NPairLiftCoupling));
      break;
    case KWInterAll:
      /* Read InterAll.def--------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NInterAll));
      break;
    case KWOneBodyG:
      /* Read cisajs.def----------------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NCisAjt));
      break;
    case KWTwoBodyG:
      /* Read cisajscktaltdc.def--------------------------------*/
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fscanf(fp,"%s %d\n", ctmp, &(X->NCisAjtCkuAlvDC));
      break;
    default:
      printf("%s", cErrIncorrectDef);
      fclose(fp);
      return -1;
      break;
    }
    /*=======================================================================*/
    fclose(fp);
  }

  if(X->iCalcModel != Spin){
    X->Ne = X->Nup+X->Ndown;
  	if(X->NLocSpn>X->Ne){
	  printf("%s", cErrNLoc);
	  printf("NLocalSpin=%d, Ne=%d\n", X->NLocSpn, X->Ne);
	  return -1;
	}
  }
  X->Nsize   = 2*X->Ne;
  X->fidx = 0;
  return 0;
}

/** 
 * 
 * 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int ReadDefFileIdxPara(
		       struct DefineList *X
		       )
{
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_CharTmpReadDef];

  int i,idx;
  int xitmp[8];
  int iKWidx=0;
  int iboolLoc=0;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  double dvalue_re, dvalue_im;
  int icnt_diagonal=0;
  int eps_CheckImag0=-12;
  eps_CheckImag0=pow(10.0, eps_CheckImag0);
  
  for(iKWidx=KWLocSpin; iKWidx< D_iKWNumDef; iKWidx++){     
    strcpy(defname, cFileNameListFile[iKWidx]);
    if(strcmp(defname,"")==0) continue;   
    printf("Read File '%s'.\n", defname);
    fp = fopen(defname, "r");
    if(fp==NULL) return ReadDefFileError(defname);
    for(i=0;i<IgnoreLinesInDef;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
    idx=0;
    /*=======================================================================*/
    switch(iKWidx){
    case KWLocSpin:
      /* Read locspn.def----------------------------------------*/
      while( fscanf(fp, "%d %d\n", &(xitmp[0]), &(xitmp[1]) )!=EOF){
	X->LocSpn[xitmp[0]] = xitmp[1];
	if(CheckSite(xitmp[1], X->Nsite) !=0){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
	idx++;
      }
      printf("Nsite= %d.\n", X->Nsite);
      if(idx!=X->Nsite){
	fclose(fp);
	return ReadDefFileError(defname);
      }
      break;
      
    case KWTrans:
      /* transfer.def--------------------------------------*/
      if(X->NTransfer>0){
	printf("X->NTransfer =%d, X->Nsite= %d.\n", X->NTransfer, X->Nsite);
	while( fscanf(fp, "%d %d %d %d %lf %lf\n", 
		      &isite1,
		      &isigma1,
		      &isite2,
		      &isigma2,
		      &dvalue_re,
		      &dvalue_im
		      )!=EOF
	       )
	  {
	    X->GeneralTransfer[idx][0]=isite1;
	    X->GeneralTransfer[idx][1]=isigma1;
	    X->GeneralTransfer[idx][2]=isite2;
	    X->GeneralTransfer[idx][3]=isigma2;
	    X->ParaGeneralTransfer[idx]=dvalue_re+dvalue_im*I;
	    if(CheckPairSite(X->GeneralTransfer[idx][0], X->GeneralTransfer[idx][2],X->Nsite) !=0){
	      fclose(fp);
	      return ReadDefFileError(defname);
	    }
	    if(isite1==isite2 && isigma1==isigma2){
	      if(abs(dvalue_im)> eps_CheckImag0){
		//NonHermite
		printf (cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, dvalue_re, dvalue_im);
		fclose(fp);
		return ReadDefFileError(defname);
	      }
	    }
	    	   
	    if(X->iCalcModel==Spin){
	      if(isite1 != isite2){
		iboolLoc=1;
		printf(cWarningIncorrectFormatForSpin2, isite1, isite2);
	      }
	    }	    
	    else if(X->iCalcModel==Kondo){
	      if(X->LocSpn[isite1]==0 || X->LocSpn[isite2]==0){
		if(isite1 != isite2){
		  iboolLoc=1;
		  printf(cErrIncorrectFormatForKondoTrans, isite1, isite2);
		}
	      }
	    }
	    idx++;
	  }
	if(idx!=X->NTransfer){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
	if(iboolLoc ==1){
	  fclose(fp);
	  return -1;
	}
      }

      if(CheckTransferHermite(X) !=0){
	printf ("%s", cErrNonHermiteTransForAll);
	fclose(fp);
	return -1;
      }
      
      break;
      
    case KWCoulombIntra:
      /*coulombintra.def----------------------------------*/
      if(X->NCoulombIntra>0){
	while( fscanf(fp, "%d %lf\n", 
		      &(X->CoulombIntra[idx][0]),
		      &(X->ParaCoulombIntra[idx])
		      )!=EOF
	       ){
	  if(CheckSite(X->CoulombIntra[idx][0], X->Nsite) !=0){
	    fclose(fp);
	    return ReadDefFileError(defname);
	  }
	  idx++;
	}
	if(idx!=X->NCoulombIntra){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;

    case KWCoulombInter:
      /*coulombinter.def----------------------------------*/
      if(X->NCoulombInter>0){
	while( fscanf(fp, "%d %d %lf\n", 
		      &(X->CoulombInter[idx][0]),
		      &(X->CoulombInter[idx][1]),
		      &(X->ParaCoulombInter[idx])
		      )!=EOF
	       ){
	  if(CheckPairSite(X->CoulombInter[idx][0], X->CoulombInter[idx][1],X->Nsite) !=0){
	    fclose(fp);
	    return ReadDefFileError(defname);
	  }

	  idx++;
	}
	if(idx!=X->NCoulombInter){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;

    case KWHund:
      /*hund.def------------------------------------------*/
      if(X->NHundCoupling>0){
	while( fscanf(fp, "%d %d %lf\n", 
		      &(X->HundCoupling[idx][0]),
		      &(X->HundCoupling[idx][1]),
		      &(X->ParaHundCoupling[idx])
		      )!=EOF
	       )
	  {
	  if(CheckPairSite(X->HundCoupling[idx][0], X->HundCoupling[idx][1],X->Nsite) !=0){
	    fclose(fp);
	    return ReadDefFileError(defname);
	  }

	  idx++;
	}
	if(idx!=X->NHundCoupling){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;
    case KWPairHop:
      /*pairhop.def---------------------------------------*/
      if(X->NPairHopping>0){
	while( fscanf(fp, "%d %d %lf\n", 
		      &(X->PairHopping[idx][0]),
		      &(X->PairHopping[idx][1]),
		      &(X->ParaPairHopping[idx])
		      )!=EOF
	       ){
	  if(CheckPairSite(X->PairHopping[idx][0], X->PairHopping[idx][1],X->Nsite) !=0){
	    fclose(fp);
	    return ReadDefFileError(defname);
	  }
	  idx++;
	}
	if(idx!=X->NPairHopping){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;

    case KWExchange:
      /*exchange.def--------------------------------------*/
      if(X->NExchangeCoupling>0){
	while( fscanf(fp, "%d %d %lf\n", 
		      &(X->ExchangeCoupling[idx][0]),
		      &(X->ExchangeCoupling[idx][1]),
		      &(X->ParaExchangeCoupling[idx])
		      )!=EOF
	       ){
	  if(CheckPairSite(X->ExchangeCoupling[idx][0], X->ExchangeCoupling[idx][1],X->Nsite) !=0){
	    fclose(fp);
	    return ReadDefFileError(defname);
	  }

	  idx++;
	}
	if(idx!=X->NExchangeCoupling){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;

    case KWIsing:
      /*ising.def--------------------------------------*/
      if(X->NIsingCoupling>0){
	while( fscanf(fp, "%d %d %lf\n", 
		      &isite1,
		      &isite2,
		      &dvalue_re
		      )!=EOF
	       ){
	  if(CheckPairSite(isite1,isite2,X->Nsite) !=0){
	    fclose(fp);
	    return ReadDefFileError(defname);
	  }

	  //input into exchange couplings
	  X->ExchangeCoupling[X->NExchangeCoupling+idx][0]=isite1;
	  X->ExchangeCoupling[X->NExchangeCoupling+idx][1]=isite2;
	  X->ParaExchangeCoupling[X->NExchangeCoupling+idx]=dvalue_re/2.0;
	  //input into inter Coulomb
	  X->CoulombInter[X->NCoulombInter+idx][0]=isite1;
	  X->CoulombInter[X->NCoulombInter+idx][1]=isite2;
	  X->ParaCoulombInter[X->NCoulombInter+idx]=-dvalue_re/4.0;
	  idx++;
	}
	if(idx!=X->NIsingCoupling){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;
      
    case KWPairLift:
      /*pairlift.def--------------------------------------*/
      if(X->NPairLiftCoupling>0){
	if(X->iCalcModel != Spin){
	  fclose(fp);
	  return -1;
	}
	while( fscanf(fp, "%d %d %lf\n", 
		      &(X->PairLiftCoupling[idx][0]),
		      &(X->PairLiftCoupling[idx][1]),
		      &(X->ParaPairLiftCoupling[idx])
		      )!=EOF
	       )
	  {
	    if(CheckPairSite(X->PairLiftCoupling[idx][0], X->PairLiftCoupling[idx][1],X->Nsite) !=0){
	      fclose(fp);
	      return ReadDefFileError(defname);
	    }

	  idx++;
	}
	if(idx!=X->NPairLiftCoupling){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;
      
    case KWInterAll:
      /*interall.def---------------------------------------*/
      X->NInterAll_Diagonal=0;
      X->NInterAll_OffDiagonal=0;
      if(X->NInterAll>0){
	while( fscanf(fp, "%d %d %d %d %d %d %d %d %lf %lf\n", 
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
		      )!=EOF
	       )
	  {
	   	    
	    if(X->iCalcModel == Spin){
	      if(!CheckFormatForSpinInt(isite1, isite2, isite3, isite4)==0){
		fclose(fp);
		return -1;
	      }
	    }
	    X->InterAll[idx][0]=isite1;
	    X->InterAll[idx][1]=isigma1;
	    X->InterAll[idx][2]=isite2;
	    X->InterAll[idx][3]=isigma2;
	    X->InterAll[idx][4]=isite3;
	    X->InterAll[idx][5]=isigma3;
	    X->InterAll[idx][6]=isite4;
	    X->InterAll[idx][7]=isigma4;
	    X->ParaInterAll[idx]=dvalue_re+I*dvalue_im;

	    if(X->InterAll[idx][0] == X->InterAll[idx][2] &&X->InterAll[idx][4] == X->InterAll[idx][6]){
	      if( X->InterAll[idx][1] == X->InterAll[idx][3] &&X->InterAll[idx][5] == X->InterAll[idx][7]){
		icnt_diagonal++;
	      }
	    }
	    if(CheckQuadSite(X->InterAll[idx][0], X->InterAll[idx][2], X->InterAll[idx][4], X->InterAll[idx][6], X->Nsite) !=0){
	      fclose(fp);
	      return ReadDefFileError(defname);
	    }
	    idx++;
	  }
	if(idx!=X->NInterAll){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }

      if(X->iCalcModel==Kondo){
	if(CheckFormatForKondoInt(X) !=0){
	  fclose(fp);
	  return -1;
	}
      }
      
      X->NInterAll_Diagonal=icnt_diagonal;
      X->NInterAll_OffDiagonal = X->NInterAll-X->NInterAll_Diagonal;
      if(!GetDiagonalInterAll(X)==0){
	fclose(fp);
	return -1;
      }

      
      if(CheckInterAllHermite(X)!=0){
	printf("%s", cErrNonHermiteInterAllForAll);
	fclose(fp);
	return -1;
      }
      
      
      break;
      
    case KWOneBodyG:
      /*cisajs.def----------------------------------------*/
      if(X->NCisAjt>0){
	while( fscanf(fp, "%d %d %d %d\n", 
		      &isite1,
		      &isigma1,
		      &isite2,
		      &isigma2) != EOF){

	  if(X->iCalcModel == Spin){
	    if(isite1 != isite2){
	      printf(cWarningIncorrectFormatForSpin2, isite1, isite2);
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
	if(idx!=X->NCisAjt){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;
      
    case KWTwoBodyG:
      /*cisajscktaltdc.def--------------------------------*/
      if(X->NCisAjtCkuAlvDC>0){
	while( fscanf(fp, "%d %d %d %d %d %d %d %d\n",
				  &isite1,
				  &isigma1,
				  &isite2,
				  &isigma2,
				  &isite3,
				  &isigma3,
				  &isite4,
				  &isigma4
				  )
		   != EOF
		   ){
		if(X->iCalcModel == Spin){
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
	if(idx!=X->NCisAjtCkuAlvDC){
	  fclose(fp);
	  return ReadDefFileError(defname);
	}
      }
      break;
    default:
      break;
    }
    fclose(fp);
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
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckSite(
	      const int iSite,
	      const int iMaxNum
	      )
{
  if(iSite>=iMaxNum) return -1;
  return 0;
}

/**
 * @brief Check Site Number for a pair -> (siteA, siteB).
 * @param[in] iSite1 a site number on a site A.
 * @param[in] iSite2 a site number on a site B.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @date 2015/07/28
 **/
int CheckPairSite(
	      const int iSite1,
	      const int iSite2,
	      const int iMaxNum
		  )
{
  if(CheckSite(iSite1, iMaxNum)!=0){
    return -1;
  }
  if(CheckSite(iSite2, iMaxNum)!=0){
    return -1;
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
    return -1;
  }
  if(CheckPairSite(iSite3, iSite4, iMaxNum)!=0){
    return -1;
  }
  return 0;
}

/**
 * @brief Check Hermite for Transfer integrals.
 * @param[in] *X Define List for getting transfer integrals.
 * @retval 0 Hermite.
 * @retval -1 NonHermite.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckTransferHermite
(
 const struct DefineList *X
 )
{
  int i,j;
  int isite1, isite2;
  int isigma1, isigma2;
  int itmpsite1, itmpsite2;
  int itmpsigma1, itmpsigma2;
  double  complex ddiff_trans;
  for(i=0; i<X->NTransfer; i++){
    isite1=X->GeneralTransfer[i][0];
    isigma1=X->GeneralTransfer[i][1];
    isite2=X->GeneralTransfer[i][2];
    isigma2=X->GeneralTransfer[i][3];
    for(j=i+1; j<X->NTransfer; j++){
      itmpsite1=X->GeneralTransfer[j][0];
      itmpsigma1=X->GeneralTransfer[j][1];
      itmpsite2=X->GeneralTransfer[j][2];
      itmpsigma2=X->GeneralTransfer[j][3];
      if(isite1 == itmpsite2 && isite2 == itmpsite1){
	if(isigma1 == itmpsigma2 && isigma2 == itmpsigma1){
	  ddiff_trans = X->ParaGeneralTransfer[i]-conj(X->ParaGeneralTransfer[j]);
	  if(cabs(ddiff_trans) > eps_CheckImag0 ){
	    printf (cErrNonHermiteTrans, isite1, isigma1, isite2, isigma2, creal(X->ParaGeneralTransfer[i]), cimag(X->ParaGeneralTransfer[i]));
	    printf (cErrNonHermiteTrans, itmpsite1, itmpsigma1, itmpsite2, itmpsigma2, creal(X->ParaGeneralTransfer[j]), cimag(X->ParaGeneralTransfer[j]));
	    return -1;
	  }
	}
      }      
    }
  }
  return 0;
}

/** 
 * 
 * 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckInterAllHermite
(
 const struct DefineList *X
 ){
  int i,j, icntincorrect, itmpret;
  int isite1, isite2, isite3, isite4;
  int isigma1, isigma2, isigma3, isigma4;
  int itmpsite1, itmpsite2, itmpsite3, itmpsite4;
  int itmpsigma1, itmpsigma2, itmpsigma3, itmpsigma4;
  double  complex ddiff_intall;
  icntincorrect=0;
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
    
    for(j=0; j<X->NInterAll_OffDiagonal; j++){
      itmpsite1=X->InterAll_OffDiagonal[j][0];
      itmpsigma1=X->InterAll_OffDiagonal[j][1];
      itmpsite2=X->InterAll_OffDiagonal[j][2];
      itmpsigma2=X->InterAll_OffDiagonal[j][3];
      itmpsite3=X->InterAll_OffDiagonal[j][4];
      itmpsigma3=X->InterAll_OffDiagonal[j][5];
      itmpsite4=X->InterAll_OffDiagonal[j][6];
      itmpsigma4=X->InterAll_OffDiagonal[j][7];
      
      if(isite1 == itmpsite4 && isite2 == itmpsite3 && isite3==itmpsite2 && isite4 == itmpsite1){
	if(isigma1 == itmpsigma4 && isigma2 == itmpsigma3 && isigma3 == itmpsigma2 && isigma4 == itmpsigma1){ 
	  itmpret=1;
	  ddiff_intall = X->ParaInterAll_OffDiagonal[i]-conj(X->ParaInterAll_OffDiagonal[j]);
	  if(cabs(ddiff_intall) > eps_CheckImag0 ){
	    printf (cErrNonHermiteInterAll, isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4, creal(X->ParaInterAll[i]), cimag(X->ParaInterAll[i]));
	    printf (cErrNonHermiteInterAll, itmpsite1, itmpsigma1, itmpsite2, itmpsigma2, itmpsite3, itmpsigma3, itmpsite4, itmpsigma4, creal(X->ParaInterAll[j]), cimag(X->ParaInterAll[j]));
	    return -1;
	  }
	}
      }
      //for spin
      else if(isite1 == itmpsite2 && isite2 ==itmpsite1 && isite3 == itmpsite4 && isite4 == itmpsite3){
	if(isigma1 == itmpsigma2 && isigma2 == itmpsigma1 && isigma3 == itmpsigma4 && isigma4 == itmpsigma3){
	  itmpret=1;
	  ddiff_intall = X->ParaInterAll_OffDiagonal[i]-conj(X->ParaInterAll_OffDiagonal[j]);
	  if(cabs(ddiff_intall) > eps_CheckImag0 ){
	    printf (cErrNonHermiteInterAll, isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4, creal(X->ParaInterAll[i]), cimag(X->ParaInterAll[i]));
	    printf (cErrNonHermiteInterAll, itmpsite1, itmpsigma1, itmpsite2, itmpsigma2, itmpsite3, itmpsigma3, itmpsite4, itmpsigma4, creal(X->ParaInterAll[j]), cimag(X->ParaInterAll[j]));
	    return -1;
	  }
	}
      }
    }
    //if counterpart for satisfying hermite conjugate does not exist.
    if(itmpret !=1){
      printf (cErrNonHermiteInterAll, isite1, isigma1, isite2, isigma2, isite3, isigma3, isite4, isigma4, creal(X->ParaInterAll_OffDiagonal[i]), cimag(X->ParaInterAll_OffDiagonal[i]));
      icntincorrect++;
    }
  }

  if( icntincorrect !=0){
    return -1;
  }  
  return 0;
}

/** 
 * 
 * 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int GetDiagonalInterAll
(
 struct DefineList *X
 )
{
  int i,icnt_diagonal, icnt_offdiagonal, tmp_i;
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
    
    if(isite1 == isite2 && isite3 == isite4){
      if(isigma1 == isigma2  && isigma3 == isigma4){
	X->InterAll_Diagonal[icnt_diagonal][0]=isite1;
	X->InterAll_Diagonal[icnt_diagonal][1]=isigma1;
	X->InterAll_Diagonal[icnt_diagonal][2]=isite3;
	X->InterAll_Diagonal[icnt_diagonal][3]=isigma3;
	X->ParaInterAll_Diagonal[icnt_diagonal] = creal(X->ParaInterAll[i]);
	icnt_diagonal++;
	continue;
      }
    }

    switch(X->iCalcModel){
    case Hubbard:
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
	  printf(cErrNonConservedInterAll,
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
	for(tmp_i=0; tmp_i<8; tmp_i++){
	  X->InterAll_OffDiagonal[icnt_offdiagonal][tmp_i]=X->InterAll[i][tmp_i];
	}
	X->ParaInterAll_OffDiagonal[icnt_offdiagonal] = X->ParaInterAll[i];
      break;
    default:
      return -1;
      break;
    }

    if(iret != -1){
      icnt_offdiagonal++;
    }
  }

  if(iret !=0){
    return -1;
  }
  
  return 0;
}

/** 
 * 
 * 
 * @param argc 
 * @param argv 
 * @param mode 
 * 
 * @return 
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
  
  if(argc == 3 && 
    (strcmp(argv[1], "-e") == 0 || 
    strcmp(argv[1], "-E") == 0  || 
    strcmp(argv[1], "--expert") == 0||
    strcmp(argv[1], "--Expert") == 0)){
    *mode=EXPERT_MODE;
  }
  else if (argc == 3 && 
    (strcmp(argv[1], "-s") ==0 || 
    strcmp(argv[1], "-S") == 0 ||
    strcmp(argv[1], "--standard") == 0 || 
    strcmp(argv[1], "--Standard") == 0 )){
    *mode=STANDARD_MODE;
  }
  else if (argc == 3 && 
    (strcmp(argv[1], "-sdry") == 0 || strcmp(argv[1], "-sDry") == 0 || 
    strcmp(argv[1], "-Sdry") == 0 || strcmp(argv[1], "-SDRY") == 0 || 
    strcmp(argv[1], "-s-dry") == 0 || strcmp(argv[1], "-S-DRY") == 0 ||
    strcmp(argv[1], "-s-Dry") == 0 || strcmp(argv[1], "-S-dry") == 0)
    ){
    *mode = STANDARD_DRY_MODE;
  }
  else{
    /*printf(cErrArgv, argv[1]);*/
    printf("\n[Usage] \n");
    printf("* Expart mode \n");
    printf("   $ HPhi -e {namelist_file} \n");
    printf("* Standard mode \n");
    printf("   $ HPhi -s {input_file} \n");
    printf("* Standard DRY mode \n");
    printf("   $ HPhi -sdry {input_file} \n");
    printf("* In this mode, Hphi stops after it generats expart input files. \n\n");
    return (-1);
  }

  return 0;
}

/** 
 * 
 * 
 * @param site1 
 * @param site2 
 * @param site3 
 * @param site4 
 * 
 * @return 
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
  else{
    printf(cWarningIncorrectFormatForSpin, site1, site2, site3, site4);
    return -1;
  }
}

/** 
 * 
 * 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckFormatForKondoInt
(
 struct DefineList *X
 ){
  int i,iboolLoc;
  int isite1, isite2, isite3, isite4;

  iboolLoc=0;
  for(i=0; i<X->NInterAll; i++){
    isite1=X->InterAll[i][0];
    isite2=X->InterAll[i][2];
    isite3=X->InterAll[i][4];
    isite4=X->InterAll[i][6];

    if(X->LocSpn[isite1]==0 || X->LocSpn[isite2]==0){
      if(isite1 != isite2){
	iboolLoc=1;
	printf(cErrIncorrectFormatForKondoInt, isite1, isite2, isite3, isite4);
	continue;
      }
    }
    if(X->LocSpn[isite3]==0 || X->LocSpn[isite4]==0){
      if(isite3 != isite4){
	iboolLoc=1;
	printf(cErrIncorrectFormatForKondoInt, isite1, isite2, isite3, isite4);
	continue;
      }
    }
  }
  if(iboolLoc==1){
    return -1;
  }  
  return 0;
}

/** 
 * 
 * 
 * @param X 
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

void ResetInteractionNum
(
 struct DefineList *X
)
{
  X->NExchangeCoupling += X->NIsingCoupling;
  X->NCoulombInter += X->NIsingCoupling;
}

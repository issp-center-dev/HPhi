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
#pragma once
#include "xsetmem.h"
#define D_FileNameMaxReadDef 256 /*!<  Max length of words for file name*/
#define D_CharTmpReadDef     200 /*!<  Max length of reading words from input files*/
#define D_CharKWDMAX     200 /*!<  Max length of words for keyword*/

/*!< Number of ignore lines in def files */
#define IgnoreLinesInDef 5

/*!< Mode of define files.*/
#define EXPERT_MODE 0
#define STANDARD_MODE 1
#define STANDARD_DRY_MODE 2

/**
 * Number of Keyword List in NameListFile for this prrogram.  
 **/
#define KWCalcMod 0
#define KWModPara 1
#define KWLocSpin 2
#define KWTrans 3
#define KWCoulombIntra 4
#define KWCoulombInter 5
#define KWHund 6
#define KWPairHop 7
#define KWExchange 8
#define KWInterAll 9
#define KWOneBodyG 10
#define KWTwoBodyG 11
#define KWPairLift 12
#define KWIsing 13
#define KWBoost 14
#define KWSingleExcitation 15
#define KWPairExcitation 16
#define KWSpectrumVec 17
#define KWLaser 18
#define KWTEOneBody 19
#define KWTETwoBody 20
#define KWThreeBodyG 21
#define KWFourBodyG 22
#define KWSixBodyG 23

int CheckSite(
          const int iListToSite,
          const int iMaxNum
          );


int CheckPairSite(
          const int iList1ToSite,
          const int iList2ToSite,
          const int iMaxNum
          );


int CheckQuadSite(
          const int iList1ToSite,
          const int iList2ToSite,
          const int iList3ToSite,
          const int iList4ToSite,
          const int iMaxNum
          );

int CheckTransferHermite
(
        struct DefineList *X
);

int CheckInterAllHermite
(
        int **InterAll,
        double complex* ParaInterAll,
        int **InterAllOffDiagonal,
        double complex*ParaInterAllOffDiagonal,
        const int NInterAllOffDiagonal,
        const int iCalcModel
);

/*
int GetDiagonalInterAll
(
 struct DefineList *X
 );
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
                                unsigned int *NChemi,
                                const int iCalcModel
                );

int JudgeDefType
(
 const int argc,
 char *argv[],
 int *mode
 );

int CheckFormatForSpinInt
(
 const int site1,
 const int site2,
 const int site3,
 const int site4
 );

/*
int CheckFormatForKondoInt
(
 struct DefineList *X
 );
*/
int CheckFormatForKondoInt
(
 const int isite1, const int isite2,
 const int isite3, const int isite4,
 int* iLocInfo
 );

int CheckFormatForKondoTrans
(
 struct DefineList *X
 );

void SetConvergenceFactor
(
 struct DefineList *X
 );

int CheckLocSpin
(
  struct DefineList *X
);

void ResetInteractionNum
(
 struct DefineList *X
 );

void InitializeInteractionNum
(
 struct DefineList *X
 );

int CheckGeneralSpinIndexForInterAll
(
 const int isite1, const int isigma1,
 const int isite2, const int isigma2,
 const int isite3, const int isigma3,
 const int isite4, const int isigma4,
 int* iLocInfo
 );

int CheckSpinIndexForTrans
(
  struct DefineList *X
 );

int CheckTotal2Sz
(
  struct DefineList *X
 );

int ReadDefFileNInt(
                    char *xNameListFile, 
                    struct DefineList *X,
                    struct BoostList *xBoost
                    );

int ReadDefFileIdxPara(
                       struct DefineList *X,
                       struct BoostList *xBoost
                       );

int CheckWords(
               const char* ctmp,
               const char* cKeyWord
               );

int GetFileNameByKW(
                    int KWidx,
                    char **FileName
                    );

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

#include "Common.h"
#include "log.h"

/**
 * @brief Functions for writing a time log.
 * @param[in] *X Define List for parameters.
 * @param[in] *cFileName character for a file name of writing.
 * @param[in] *cTimeKeeper_Message character for a sentence writing in a file.
 * @param[in] *cWriteType character for determining writing option for fopen.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int TimeKeeper
(
 struct BindStruct *X,
 const char *cFileName,
 const char *cTimeKeeper_Message,
 const char *cWriteType
 )
{
  FILE *fp;
  char sdt[D_FileNameMax];
  struct tm *area;
  time_t tx;
  
  sprintf(sdt, cFileName, X->Def.CDataFileHead);
  tx   = time(NULL);
  area = localtime(&tx);
  if(childfopenMPI(sdt, cWriteType, &fp)!=0){
    return -1;
  }
  fprintf(fp, cTimeKeeper_Message, asctime(area));
  fcloseMPI(fp);
  return 0;
}

/**
 * @brief Functions for writing a time log.
 * @param[in] *X Define List for parameters.
 * @param[in] *cFileName character for a file name of writing.
 * @param[in] *cTimeKeeper_Message character for a sentence writing in a file.
 * @param[in] *cWriteType character for determining writing option for fopen.
 * @param[in] istep int for writing steps of progress for a function under calculating.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int TimeKeeperWithStep
(
 struct BindStruct *X,
 const char *cFileName,
 const char *cTimeKeeper_Message,
 const char *cWriteType,
 const int istep
 )
{
  FILE *fp;
  char sdt[D_FileNameMax];
  struct tm *area;
  time_t tx;
  
  sprintf(sdt, cFileName, X->Def.CDataFileHead);
  tx   = time(NULL);
  area = localtime(&tx);
  if(childfopenMPI(sdt, cWriteType, &fp)!=0){
    return -1;
  }
  fprintf(fp, cTimeKeeper_Message, istep, asctime(area));
  fcloseMPI(fp);
  return 0;
}


/**
 * @brief Functions for writing a time log.
 * @param[in] *X Define List for parameters.
 * @param[in] *cFileName character for a file name of writing.
 * @param[in] *cTimeKeeper_Message character for a sentence writing in a file.
 * @param[in] *cWriteType character for determining writing option for fopen.
 * @param[in] irand int for writing a number of seed of rand.
 * @param[in] istep int for writing steps of progress for a function under calculating.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int TimeKeeperWithRandAndStep
(
 struct BindStruct *X,
 const char *cFileName,
 const char *cTimeKeeper_Message,
 const char *cWriteType,
 const int irand,
 const int istep
 ){
  FILE *fp;
  char sdt[D_FileNameMax];
  struct tm *area;
  time_t tx;
  
  sprintf(sdt, cFileName, X->Def.CDataFileHead);
  tx   = time(NULL);
  area = localtime(&tx);
  if(childfopenMPI(sdt, cWriteType, &fp)!=0){
    return -1;
  }
  fprintf(fp, cTimeKeeper_Message, irand, istep, asctime(area));
  fcloseMPI(fp);
  return 0;
  
}

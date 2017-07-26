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
/**@file
@brief Functions to open file(s) in output/ directory.
*/
#include "FileIO.h"
#include "wrapperMPI.h"
/** 
@brief Only the root process open file in output/ directory
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@return -1 if file cannot be opened, 0 if not.
*/
int childfopenMPI(
  const char* _cPathChild,//!<[in] File name
  const char* _cmode,//!<[in] "w", "r", etc.
  FILE **_fp//!<[inout] File pointer
){
  char ctmpPath[D_FileNameMax]="";
  strcat(ctmpPath, cParentOutputFolder);
  strcat(ctmpPath, _cPathChild);
  *_fp = fopenMPI(ctmpPath, _cmode);
  
  if(*_fp == NULL){
    fprintf(stdoutMPI, cErrFIOpen, ctmpPath);
    return -1;
  }
    
  return 0;
}/*int childfopenMPI*/
/** 
@brief All processes open file in output/ directory
@author Takahiro Misawa (The University of Tokyo)
@author Kazuyoshi Yoshimi (The University of Tokyo)
@return -1 if file cannot be opened, 0 if not.
*/
int childfopenALL(
  const char* _cPathChild,//!<[in] File name
  const char* _cmode,//!<[in] "w", "r", etc.
  FILE **_fp//!<[inout] File pointr
) {
  char ctmpPath[D_FileNameMax]="";
  strcat(ctmpPath, cParentOutputFolder);
  strcat(ctmpPath, _cPathChild);
  *_fp = fopen(ctmpPath, _cmode);
  
  if(*_fp == NULL){
    fprintf(stdout, cErrFIOpen, ctmpPath);
    return -1;
  }
    
  return 0;
}/*int childfopenALL*/

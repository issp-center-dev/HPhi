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
#include "FileIO.h"
#include "wrapperMPI.h"

/** 
 * @brief 
 * 
 * @param[in] _cPathChild 
 * @param[in] _cmode 
 * @param[in] _fp 
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @return 
 */
int childfopenMPI(const char* _cPathChild, const char* _cmode, FILE **_fp){
  char ctmpPath[D_FileNameMax]="";
  strcat(ctmpPath, cParentOutputFolder);
  strcat(ctmpPath, _cPathChild);
  *_fp = fopenMPI(ctmpPath, _cmode);
  
  if(*_fp == NULL){
    fprintf(stdoutMPI, cErrFIOpen, ctmpPath);
    return -1;
  }
    
  return 0;
}

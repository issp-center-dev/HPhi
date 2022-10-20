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
#include "global.h"
#include "LogMessage.h"
#include "struct.h"

int TimeKeeper
(
 struct BindStruct *X,
 const char *cFileName,
 const char *cTimeKeeper_Message,
 const char *cWriteType
);

int TimeKeeperWithStep
(
 struct BindStruct *X,
 const char *cFileName,
 const char *cTimeKeeper_Message,
 const char *cWriteType,
 const int istep
 );


int TimeKeeperWithRandAndStep
(
 struct BindStruct *X,
 const char *cFileName,
 const char *cTimeKeeper_Message,
 const char *cWriteType,
 const int istep,
 const int irand
 );

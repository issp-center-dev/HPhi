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
/*-------------------------------------------------------------*/
/**@file
@brief Print logo mark and version number
*/
#include <stdio.h>
#include "global.h"
/**
@brief Print logo mark and version number
*/
void splash(){
   int ver_maj =
#include "version_major.h"
;
   int ver_min =
#include "version_minor.h"
;
   int ver_pat =
#include "version_patch.h"
;

  fprintf(stdoutMPI, "                                                                \n");
  fprintf(stdoutMPI, "      ,ammmmmmmmmmmmmmb,,        Welcome to the                 \n");
  fprintf(stdoutMPI, "    ,@@` dm          mb  ===m                                    \n");
  fprintf(stdoutMPI, "  ,@@` d@@@@@@@@@@@@@@@@b Pm,    @@          @@        @@       \n");
  fprintf(stdoutMPI, " d@  d@@@ @@@ @@@@@@ @@@@b ~@a   @@          @@     @@@@@@@@    \n");
  fprintf(stdoutMPI, "d@   @@@@ ^^^ @@@@ m m @@@   @,  @@          @@   @@@  @@  @@@  \n");
  fprintf(stdoutMPI, "@    @@@@_@@@_@@@@mm mm@@@   @|  @@mmmmmmmmmm@@  @@    @@    @@ \n");
  fprintf(stdoutMPI, "P@    9@@@@@@@@@@@@@@@@@P    @~  @@@@@@@@@@@@@@  @@    @@    @@ \n");
  fprintf(stdoutMPI, " @@      ~~9@@@@@@PPP~      @P   @@          @@   @@@  @@  @@@  \n");
  fprintf(stdoutMPI, "  ~@@b      @@@@@@@      ,@@~    @@          @@     @@@@@@@@    \n");
  fprintf(stdoutMPI, "    ~@@@m,,@@@@@@@@@  ,m@~`      @@          @@        @@       \n");
  fprintf(stdoutMPI, "        ~~9@@@@@@@@@  ~                                         \n");
  fprintf(stdoutMPI, "           9@P~~~9@P             Version %d.%d.%d    \n", ver_maj, ver_min, ver_pat);
  fprintf(stdoutMPI, "                                                                \n");

}/*void splash()*/

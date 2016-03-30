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
/*-------------------------------------------------------------*/

#include <stdio.h>
#include "wrapperMPI.h"

void splash(){

  fprintf(stdoutMPI, "                                                                \n");
  fprintf(stdoutMPI, "      ,ammmmmmmmmmmmmmb,,        Welcome to the                 \n");
  fprintf(stdoutMPI, "    ,@@` mmmmmmmmmmmmmm ===m                                    \n");
  fprintf(stdoutMPI, "  ,@@` d@@@@@@@@@@@@@@@@b Pm,    @@          @@        @@       \n");
  fprintf(stdoutMPI, " d@  d@@@ @@@ @@@@@@ @@@@b ~@a   @@          @@     @@@@@@@@    \n");
  fprintf(stdoutMPI, "d@   @@@@ ^^^ @@@@ m m @@@   @,  @@          @@   @@@  @@  @@@  \n");
  fprintf(stdoutMPI, "@    @@@@_@@@_@@@@mm mm@@@   @|  @@mmmmmmmmmm@@  @@    @@    @@ \n");
  fprintf(stdoutMPI, "P@    9@@@@@@@@@@@@@@@@@P    @~  @@@@@@@@@@@@@@  @@    @@    @@ \n");
  fprintf(stdoutMPI, " @@      ~~9@@@@@@PPP~      @P   @@          @@   @@@  @@  @@@  \n");
  fprintf(stdoutMPI, "  ~@@b      @@@@@@@      ,@@~    @@          @@     @@@@@@@@    \n");
  fprintf(stdoutMPI, "    ~@@@m,,@@@@@@@@@  ,m@~`      @@          @@        @@       \n");
  fprintf(stdoutMPI, "        ~~9@@@@@@@@@  ~                                         \n");
  fprintf(stdoutMPI, "           9@P~~~9@P             Version 1.0                    \n");
  fprintf(stdoutMPI, "                                                                \n");

}

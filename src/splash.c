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

  fprintf(stdoutMPI, "                                                                                                              \n");
  fprintf(stdoutMPI, "                  ,,m@@@@@@@@@@@@@@@mm..                                                                      \n");
  fprintf(stdoutMPI, "              ,m@@@@PPP~~~~~~~~~~~~PP@@@@@m,,               Welcome to the                                    \n");
  fprintf(stdoutMPI, "          ,m@@@@ ,,,  ,mm@@@@@@@@mm,,   ~~@@@@m                                                               \n");
  fprintf(stdoutMPI, "        ,@@@P  /@@@@@@@@@@@@@@@@@@@@@@@@mmmm, @@@m                                                            \n");
  fprintf(stdoutMPI, "      ,@@@    |@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   @@@@                                                          \n");
  fprintf(stdoutMPI, "     @@@    ,@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.    @@@       @@@              @@@                @@@           \n");
  fprintf(stdoutMPI, "    @@P    m@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     @@,     @@@              @@@                @@@           \n");
  fprintf(stdoutMPI, "   @@     f@@@@@@@@@@@@@@@@@@@@@@@@@@P9@@@@@@@@@     @@,    @@@              @@@          @@@@@@@@@@@@@@@     \n");
  fprintf(stdoutMPI, "  @@     @@@@@@@@~~@@@@@P9@@@@@@@@P~~  ~~@@@@@@@@     @@.   @@@              @@@        @@@     @@@     @@@   \n");
  fprintf(stdoutMPI, " @@P     @@@@@@@m  @@@@@  @@@@@@P  ,m  m, ~@@@@@@      @@   @@@              @@@       @@@      @@@      @@@  \n");
  fprintf(stdoutMPI, " @@      @@@@@@@@         @@@@@@|  @@  @@ ,@@@@@@      @@   @@@@@@@@@@@@@@@@@@@@      @@@       @@@       @@@ \n");
  fprintf(stdoutMPI, " @@      @@@@@@@~ ,@@@@@  @@@@@@@,       ,@@@@@@P      @@   @@@@@@@@@@@@@@@@@@@@      @@@       @@@       @@@ \n");
  fprintf(stdoutMPI, " @@       @@@@@@@m@@@@@@mm@@@@@@@@@@b  m@@@@@@@P       @@   @@@              @@@       @@@      @@@      @@@  \n");
  fprintf(stdoutMPI, " @@L       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@P        @@   @@@              @@@        @@@     @@@     @@@   \n");
  fprintf(stdoutMPI, "  @@         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@~         @@P   @@@              @@@          @@@@@@@@@@@@@@@     \n");
  fprintf(stdoutMPI, "  ~@@           @@@@@@@@@@@@@@@@@@@@@@@@P~           @@P    @@@              @@@                @@@           \n");
  fprintf(stdoutMPI, "   ~@@               @@@@@@@@@@@@@P~~               @@@     @@@              @@@                @@@           \n");
  fprintf(stdoutMPI, "    ~@@,                @@@@@@@@@@L               ,@@P                                                        \n");
  fprintf(stdoutMPI, "     ~@@@,             @@@@@@@@@@@@@            ,@@@~       Version 0.3                                       \n");
  fprintf(stdoutMPI, "       ~@@@m,         @@@@@@@@@@@@@@L         ,@@@~                                                           \n");
  fprintf(stdoutMPI, "         ~9@@@m,     @@@@@@@@@@@@@@@@     .@@@@@~                                                             \n");
  fprintf(stdoutMPI, "            ~P@@@@mm,@@@@@@@@@@@@@@@@| @@@@@@~                                                                \n");
  fprintf(stdoutMPI, "                ~@@@@@@@@@@@@@@@@@@@@@ ~P~                                                                    \n");
  fprintf(stdoutMPI, "                     |@@@@@@@@@@@@@@@                                                                         \n");
  fprintf(stdoutMPI, "                      @@@@@P^^^^9@@@P                                                                         \n");
  fprintf(stdoutMPI, "                       ~~~       ~~~                                                                          \n");
  fprintf(stdoutMPI, "                                                                                                              \n");

}

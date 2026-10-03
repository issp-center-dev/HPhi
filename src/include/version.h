/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef HPHI_VERSION_H
#define HPHI_VERSION_H

/*
The version number is defined only here.
CMakeLists.txt, dist.sh and doc/(en|ja)/source/conf.py read the three lines
below, so keep the form "#define NAME number".
*/
#define HPHI_VERSION_MAJOR 3
#define HPHI_VERSION_MINOR 7
#define HPHI_VERSION_PATCH 0

const char *GetGitHash(void);

#endif

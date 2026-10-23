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
#include "version.h"

/*
version_git.h is written into the build directory by cmake/git_hash.cmake.
*/
#ifdef HPHI_HAVE_VERSION_GIT_H
#include "version_git.h"
#endif
#ifndef HPHI_GIT_HASH
#define HPHI_GIT_HASH ""
#endif

/**
@brief Abbreviated hash (8 digits) of the commit which HPhi was built from
@return The hash, followed by "-dirty" if the source had changes which were
not committed. An empty string if the commit is not known.
*/
const char *GetGitHash(void) {
  return HPHI_GIT_HASH;
}/*const char *GetGitHash*/

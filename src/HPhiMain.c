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

#include <sz.h>
#include <HPhiTrans.h>
#include <output_list.h>
#include <diagonalcalc.h>
#include <CalcByLanczos.h>
#include <CalcByLOBPCG.h>
#include <CalcByFullDiag.h>
#include <CalcByTPQ.h>
#include <CalcSpectrum.h>
#include <check.h>
#include "CalcByCanonicalTPQ.h"
#include "CalcByTEM.h"
#include "readdef.h"
#include "StdFace_main.h"
#include "wrapperMPI.h"
#include "splash.h"
#include "CalcTime.h"

/*!
  @mainpage

  <H2>Introduction</H2>
  A numerical solver package for a wide range of quantum lattice models including Hubbard-type itinerant electron hamiltonians, quantum spin models, and Kondo-type hamiltonians for itinerant electrons coupled with quantum spins. The Lanczos algorithm for finding ground states and newly developed Lanczos-based algorithm for finite-temperature properties of these models are implemented for parallel computing. A broad spectrum of users including experimental researchers is cordially welcome.
  <HR>
  <H2>Developers</H2>
  Youhei Yamaji (Quantum-Phase Electronics Center, The University of Tokyo)\n
  Takahiro Misawa (Institute for Solid State Physics, The University of Tokyo)\n
  Synge Todo (Department of Physics, The University of Tokyo)\n
  Kazuyoshi Yoshimi (Institute for Solid State Physics, The University of Tokyo)\n
  Mitsuaki Kawamura (Institute for Solid State Physics, The University of Tokyo)\n
  Kota Ido (Department of Applied Physics, The University of Tokyo)\n
  Naoki Kawashima (Institute for Solid State Physics, The University of Tokyo)
  <HR>
  <H2>Methods</H2>
  - Lanczos algorithm
  - Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG) method : See LOBPCG_Main()
  - Thermal Pure Quantum (TPQ) state
  - Full Diagonalization
  - Shifted BiCG method : See CalcSpectrumByBiCG()
  - Lehmann's spectral representation : See CalcSpectrumByFullDiag()
  .
  <HR>
  <H2>Target models</H2>
  Hubbard model, Heisenberg model, Kondo lattice model, Kitaev model, Kitaev-Heisenberg model, multi-orbital Hubbard model
  <HR>
  <H2>Important functions and source files</H2>
  <ul>
    <li>mltply.c : Perform Hamiltonian-vector product</li>
    <ul>
      <li>mltplyHubbard.c : For Hubbard and Kondo system</li>
      <li>mltplySpin.c : For local spin system</li>
    </ul>
    <li>StdFace_main.c : Construct typical models</li>
    <li>global.h : Global variables</li>
    <li>struct.h : Binded struct</li>
  </ul>
  <HR>
  <H2>How to modify HPhi (Developer's note)</H2>
  - @ref page_codingrule
  - @ref page_addstandard
  - @ref page_addstandardval
  - @ref page_variable
  - @ref page_setmem
  - @ref page_cmake
  - @ref page_addcalcmod
  - @ref page_addmodpara
  - @ref page_addexpert
  - @ref page_time
  - @ref page_log
  - Some contrivances for HPhi (only in Jananese) https://www.pasums.issp.u-tokyo.ac.jp/wp-content/themes/HPhi/media/develop/tips.pdf
  .
  <HR>
  <H2>Link</H2>
  https://github.com/issp-center-dev/HPhi
  <HR>
  <H2>Download</H2>
  https://github.com/issp-center-dev/HPhi/releases
  <HR>
  <H2>Forum</H2>
  https://github.com/issp-center-dev/HPhi/issues
  <HR>
  <H2>licence</H2>
  <B>GNU GPL version 3</B>\n
  This software is developed under the support of "Project for advancement of software usability in materials science" by The Institute for Solid State Physics, The University of Tokyo.\n

@page page_codingrule Coding rule

- Do not use TAB character. Use two spaces as an indent.
- Use @c default(none) for scoping of OpenMP-parallel region. E.g.
  @dontinclude CalcByLOBPCG.c
  @skip pragma
  @until 0.0
- Variable declared with @c const must not be included in @c firstprivate of OpenMP scoping.
  Use @c shared.
- For MPI parallelization, use the following functions for I/O and abortation:
  - fgetsMPI() instead of @c fgets
  - @c fprintf(::stdoutMPI,... instead of @c printf(...
  - fopenMPI() instead of @c fopen
  - exitMPI() instead of @c exit
- When you add new features into HPhi, please run <tt>make test</tt>,
  and check whether other features still work fine.
  Also, try <tt>make test MPIRUN="mpiexec -np 4"</tt> to check MPI feature.
.    

@page page_cmake Add new source-file, executable, scripts (handle CMake)

To build HPhi, CMake is required.
We have to modify the CMake configuration file when we add new sources, executables, scripts.

@section sec_newsource New source code

When we add a new source code, we have to add the file-name
into the following part of @c src/CMakeLists.txt.

\code{cmake}
set(SOURCES source1.c source2.c ...)
\endcode

@section sec_newexecutable New executable

When we add a new executable ("myprog" in this case),
we have to add following command in @c src/CMakeLists.txt.

\code{CMake}
set(SOURCES_MYPROG source1.c source2.c ...)
add_executable(myprog ${SOURCES_MYPROG})
target_link_libraries(myprog ${LAPACK_LIBRARIES} m)
if(MPI_FOUND)
target_link_libraries(myprog ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
install(TARGETS myprog RUNTIME DESTINATION bin)
\endcode

@section sec_newscript New script

When we add a new script written in python, sh, etc. ("myscript.sh" in this case)
into @c tool/, we have to add the following command in @c tool/CMakeLists.txt.

\code{CMake}
configure_file(myscript.sh myscript.sh COPYONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/myscript.sh DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
\endcode
*/

/** 
 * @brief Main program for HPhi
 * 
 * @param argc [in] argument count
 * @param argv [in] argument vector
 *
 * @version 2.1 Add Time evolution mode.
 * @version 1.2 Add calculation spectrum mode.
 * @version 1.0
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @retval -1 fail the calculation.
 * @retval 0 succeed the calculation.
 */
int main(int argc, char* argv[]){

  int mode=0;
  char cFileListName[D_FileNameMax];

  stdoutMPI = stdout;
  if(JudgeDefType(argc, argv, &mode)!=0){
      exitMPI(-1);
  }

  if (mode == STANDARD_DRY_MODE) {
    myrank = 0;
    nproc = 1;
    stdoutMPI = stdout;
    splash();
  }
  else InitializeMPI(argc, argv);

  //Timer
  InitTimer();
  if (mode != STANDARD_DRY_MODE) StartTimer(0);
 
  //MakeDirectory for output
  struct stat tmpst;
  if (myrank == 0) {
    if (stat(cParentOutputFolder, &tmpst) != 0) {
      if (mkdir(cParentOutputFolder, 0777) != 0) {
        fprintf(stdoutMPI, "%s", cErrOutput);
        exitMPI(-1);
      }
    }
  }/*if (myrank == 0)*/

  strcpy(cFileListName, argv[2]);
  
  if(mode==STANDARD_MODE || mode == STANDARD_DRY_MODE){
    if (myrank == 0) StdFace_main(argv[2]);
    strcpy(cFileListName, "namelist.def");
    if (mode == STANDARD_DRY_MODE){
      fprintf(stdout, "Dry run is Finished. \n\n");
      return 0;
    }
  }

  setmem_HEAD(&X.Bind);
  if(ReadDefFileNInt(cFileListName, &(X.Bind.Def), &(X.Bind.Boost))!=0){
    fprintf(stdoutMPI, "%s", cErrDefFile);
    exitMPI(-1);
  }

  if (X.Bind.Def.nvec < X.Bind.Def.k_exct){
    fprintf(stdoutMPI, "%s", cErrnvec);
    fprintf(stdoutMPI, cErrnvecShow, X.Bind.Def.nvec, X.Bind.Def.k_exct);
    exitMPI(-1);
  }
  fprintf(stdoutMPI, "%s", cProFinishDefFiles);
  
  /*ALLOCATE-------------------------------------------*/
  setmem_def(&X.Bind, &X.Bind.Boost);
  /*-----------------------------------------------------*/

  /*Read Def files.*/
  TimeKeeper(&(X.Bind), cFileNameTimeKeep, cReadDefStart, "w");
  if(ReadDefFileIdxPara(&(X.Bind.Def), &(X.Bind.Boost))!=0){
    fprintf(stdoutMPI, "%s", cErrIndices);
    exitMPI(-1);
  }
  TimeKeeper(&(X.Bind), cFileNameTimeKeep, cReadDefFinish, "a");
  fprintf(stdoutMPI, "%s", cProFinishDefCheck);

  /*Set convergence Factor*/
  SetConvergenceFactor(&(X.Bind.Def));

  if (X.Bind.Def.iCalcType == FullDiag
      && X.Bind.Def.iFlgScaLAPACK ==0
      && nproc != 1) {
    fprintf(stdoutMPI, "Error: Full Diagonalization by LAPACK is only allowed for one process.\n");
    FinalizeMPI();
    return(-1);
  }

  /*---------------------------*/
  if(HPhiTrans(&(X.Bind))!=0) {
    exitMPI(-1);
  }

  //Start Calculation
  if(X.Bind.Def.iFlgCalcSpec == CALCSPEC_NOT) {
    
    if(check(&(X.Bind))==MPIFALSE){
     exitMPI(-1);
    }
    
    /*LARGE VECTORS ARE ALLOCATED*/
    if (setmem_large(&X.Bind) != 0) {
      fprintf(stdoutMPI, cErrLargeMem, iErrCodeMem);
      exitMPI(-1);
    }
    
    StartTimer(1000);
    if(sz(&(X.Bind), list_1, list_2_1, list_2_2)!=0){
      exitMPI(-1);
    }

    StopTimer(1000);
    if(X.Bind.Def.WRITE==1){
      output_list(&(X.Bind));
      exitMPI(-2);
    }
    StartTimer(2000);
    diagonalcalc(&(X.Bind));
    StopTimer(2000);
      
    switch (X.Bind.Def.iCalcType) {
    case Lanczos:
      StartTimer(4000);
      if (CalcByLanczos(&X) != TRUE) {
        StopTimer(4000);
        exitMPI(-3);
      }
      StopTimer(4000);
      break;

    case CG:
      if (CalcByLOBPCG(&X) != TRUE) {
          exitMPI(-3);
      }
      break;

      case FullDiag:
        StartTimer(5000);
        if (CalcByFullDiag(&X) != TRUE) {
            FinalizeMPI();
        }
        StopTimer(5000);
      break;

      case TPQCalc:
        StartTimer(3000);        
        if (CalcByTPQ(NumAve, X.Bind.Def.Param.ExpecInterval, &X) != TRUE) {
          StopTimer(3000);
          exitMPI(-3);
        }
        StopTimer(3000);
      break;

      case cTPQ:
        StartTimer(3000);        
        if (CalcByCanonicalTPQ(NumAve, X.Bind.Def.Param.ExpecInterval, &X) != TRUE) {
          StopTimer(3000);
          exitMPI(-3);
        }
        StopTimer(3000);
      break;


      case TimeEvolution:
        if(CalcByTEM(X.Bind.Def.Param.ExpecInterval, &X)!=0){
            exitMPI(-3);
        }
      break;

    default:
      StopTimer(0);
      exitMPI(-3);
    }
  }
  else{
    StartTimer(6000);
    if (CalcSpectrum(&X) != TRUE) {
      StopTimer(6000);
      exitMPI(-3);
    }
    StopTimer(6000);
  }
  
  StopTimer(0);
  OutputTimer(&(X.Bind));
  FinalizeMPI();
  return 0;
}

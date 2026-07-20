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
#include "symmetry_basis.h"
#include "symmetry_basis_io.h"
#include "symmetry_matvec_plan.h"

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
  - Add new lattice model into Standard mode (see StdFace documentation)
  - Add new input variable into Standard mode (see StdFace documentation)
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

@section sec_coding_general General Rules

- Do not use TAB character. Use two spaces as an indent.
- Use C99 standard for compatibility.
- All source files should include the GPL license header.
- Use Doxygen-style comments for all public functions.

@section sec_coding_naming Naming Conventions

@subsection subsec_functions Function Names
- Use CamelCase with descriptive prefixes:
  - @c child_ : Internal functions called by main routines
  - @c X_child_ : Internal functions with explicit struct parameter
  - @c _GetInfo : Functions that retrieve information
  - @c _MPI : MPI-specific implementations
  - @c _MPIsingle / @c _MPIdouble : Single/double MPI communication patterns

Example: @c X_child_GC_general_hopp_MPIsingle()

@subsection subsec_variables Variable Names
- Loop indices: @c i, @c j, @c idim, @c jdim
- Site indices: @c isite1, @c isite2
- Spin indices: @c sigma1, @c sigma2, @c ispin
- Bit representations: @c ibit, @c ibitsite
- Temporary values: @c tmp_v0, @c tmp_v1, @c tmp_trans
- Damping/product: @c dam_pr (Hamiltonian matrix element)

@subsection subsec_macros Macro Names
- Use ALL_CAPS with @c D_ prefix for constants:
  - @c D_FileNameMax
  - @c D_CharTmpReadDef
- Model types: @c Hubbard, @c Spin, @c HubbardGC, @c SpinGC, etc.
- Calculation types: @c Lanczos, @c TPQCalc, @c FullDiag, @c CG

@section sec_coding_openmp OpenMP Parallelization

- Always use @c default(none) for scoping of OpenMP-parallel region:
  @dontinclude CalcByLOBPCG.c
  @skip pragma
  @until 0.0
- Classify all variables explicitly:
  - @c private() : Thread-local variables (loop indices)
  - @c shared() : Read-only shared data, const variables
  - @c firstprivate() : Variables initialized from master thread
  - @c reduction() : Variables for accumulation (+, *, max, min)
- Variable declared with @c const must not be included in @c firstprivate.
  Use @c shared instead.
- Use conditional compilation for OpenMP headers:
\code{c}
#ifdef _OPENMP
#include <omp.h>
#endif
\endcode

@section sec_coding_mpi MPI Parallelization

For MPI parallelization, use the following wrapper functions:
- fgetsMPI() instead of @c fgets
- @c fprintf(::stdoutMPI,... instead of @c printf(...
- fopenMPI() instead of @c fopen
- exitMPI() instead of @c exit

Additional MPI utilities in wrapperMPI.h:
- MaxMPI_li(), MaxMPI_d() : Get maximum across processes
- SumMPI_li(), SumMPI_dc() : Sum across processes
- NormMPI_dc() : Compute vector norm
- VecProdMPI() : Compute vector inner product

Use conditional compilation for MPI-specific code:
\code{c}
#ifdef MPI
// MPI-specific implementation
#endif
\endcode

@section sec_coding_memory Memory Management

- Use dedicated allocation functions from setmemory.c:
  - @c i_1d_allocate(), @c i_2d_allocate() : int arrays
  - @c d_1d_allocate(), @c d_2d_allocate() : double arrays
  - @c cd_1d_allocate(), @c cd_2d_allocate() : complex double arrays
  - @c lui_1d_allocate() : unsigned long int arrays
- Always use @c calloc (not @c malloc) for automatic zero-initialization.
- Check allocation results for critical allocations.

@section sec_coding_error Error Handling

- Return @c int for error codes: 0 = success, -1 = error
- Use global error messages defined in ErrorMessage.h/c
- Always check file operations:
\code{c}
fp = fopenMPI(filename, "r");
if (fp == NULL) {
  fprintf(stdoutMPI, cErrFIOpen, filename);
  return -1;
}
\endcode

@section sec_coding_testing Testing

When you add new features into HPhi, please run:
\code{bash}
make test
\endcode
and check whether other features still work fine.
Also, try MPI tests:
\code{bash}
make test MPIRUN="mpiexec -np 4"
\endcode

@page page_cmake Add new source-file, executable, scripts (handle CMake)

To build HPhi, CMake is required.
We have to modify the CMake configuration file when we add new sources, executables, scripts.

@section sec_newsource New source code

When we add a new source code, we have to add the file-name
into the following part of @c src/CMakeLists.txt.

\code{cmake}
set(SOURCES source1.c source2.c ...)
\endcode

If you also add a header file, place it in @c src/include/ directory.
Header files are automatically found by CMake.

@section sec_newexecutable New executable

When we add a new executable ("myprog" in this case),
we have to add following command in @c src/CMakeLists.txt.

\code{cmake}
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

\code{cmake}
configure_file(myscript.sh myscript.sh COPYONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/myscript.sh DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
\endcode

@section sec_cmake_options CMake Build Options

HPhi supports several CMake options:
- @c -DUSE_SCALAPACK=ON : Enable ScaLAPACK for parallel full diagonalization
- @c -DCMAKE_C_COMPILER=mpicc : Specify MPI compiler
- @c -DCMAKE_BUILD_TYPE=Release : Build with optimizations

Example build commands:
\code{bash}
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
make install
\endcode

@page page_variable Global variables and Data structure

@section sec_global_vars Global Variables

Global variables are defined in global.h and global.c.

@subsection subsec_mpi_vars MPI-related Variables
- @c nproc : Total number of MPI processes
- @c myrank : Rank of current MPI process
- @c nthreads : Number of OpenMP threads
- @c stdoutMPI : MPI-safe stdout (only rank 0 writes)

@subsection subsec_path_vars Path and File Variables
- @c cParentOutputFolder : Output directory path ("output/")
- @c cFileNameTimeKeep : Timer log filename

@subsection subsec_wave_vars Wavefunction Vectors
- @c v0, @c v1 : Main wavefunction vectors
- @c v1buf : MPI communication buffer
- @c list_1 : Hilbert space basis list
- @c list_2_1, @c list_2_2 : Two-body basis lists

@section sec_data_struct Main Data Structures

@subsection subsec_bindstruct BindStruct (X)
The main structure containing all calculation data:
\code{c}
struct BindStruct {
  struct DefineList Def;   // Input parameters and definitions
  struct CheckList Check;  // Dimension checks
  struct LargeList Large;  // Large temporary data
  struct PhysList Phys;    // Physical quantities
  struct BoostList Boost;  // Boost optimization data
};
\endcode

@subsection subsec_definelist DefineList (X->Def)
Contains all input parameters:
- @c Nsite : Number of sites
- @c Ne, @c NeMPI : Number of electrons (local/total)
- @c Nup, @c Ndown : Number of up/down spins
- @c iCalcType : Calculation method (Lanczos, TPQ, etc.)
- @c iCalcModel : Model type (Hubbard, Spin, etc.)
- @c NTransfer : Number of transfer terms
- @c NInterAll : Number of interaction terms
- @c GeneralTransfer[][] : Transfer integral indices
- @c ParaGeneralTransfer[] : Transfer integral values

@subsection subsec_largelist LargeList (X->Large)
Temporary data for Hamiltonian operations:
- @c prdct : Inner product accumulator
- @c is1_spin, @c is2_spin : Spin masks
- @c isA_spin, @c isB_spin : Additional spin masks
- @c tmp_trans : Current transfer value
- @c mode : Operation mode flag

@subsection subsec_physlist PhysList (X->Phys)
Physical quantities:
- @c energy : Ground state energy
- @c var : Energy variance
- @c doublon : Double occupancy
- @c num : Particle number
- @c Sz : Total \f$S_z\f$
- @c s2 : Total \f$S^2\f$

@page page_log Various kind of Logs

@section sec_log_progress Progress Messages

Progress messages are defined in ProgressMessage.h/c.
Use these predefined strings for consistent output:
\code{c}
fprintf(stdoutMPI, "%s", cProFinishDefFiles);  // "Read definition files."
fprintf(stdoutMPI, "%s", cProFinishDefCheck);  // "Check definition files."
\endcode

@section sec_log_error Error Messages

Error messages are defined in ErrorMessage.h/c.
Examples:
- @c cErrNameList : Command line argument error
- @c cErrDefFile : Definition file error
- @c cErrOutput : Output directory creation error
- @c cErrLargeMem : Memory allocation error

@section sec_log_time Time Logging

Use TimeKeeper() for elapsed time logging:
\code{c}
TimeKeeper(&(X.Bind), cFileNameTimeKeep, cReadDefStart, "w");
// ... processing ...
TimeKeeper(&(X.Bind), cFileNameTimeKeep, cReadDefFinish, "a");
\endcode

Output is written to @c output/Time_*.dat files.

@section sec_log_debug Debug Output

For debugging, use conditional output:
\code{c}
#ifdef DEBUG
fprintf(stdoutMPI, "Debug: value = %d\n", value);
#endif
\endcode

@page page_addexpert Add new input-file for Expert mode

@section sec_expert_overview Overview

Expert mode input files are processed in readdef.c.
Each input file type is registered in the file list (namelist.def).

@section sec_expert_steps Steps to Add New Input File

@subsection subsec_expert_step1 Step 1: Define Keywords
Add keywords to DefCommon.h:
\code{c}
#define KWMyNewFile 100  // Unique keyword ID
\endcode

@subsection subsec_expert_step2 Step 2: Add File Type
In readdef.c, add to cKWListOfFileNameList[]:
\code{c}
cKWListOfFileNameList[KWMyNewFile] = "mynewfile";
\endcode

@subsection subsec_expert_step3 Step 3: Implement Reader
Create a new reader function:
\code{c}
int ReadMyNewFile(const char* filename, struct DefineList *Def) {
  FILE *fp = fopenMPI(filename, "r");
  if (fp == NULL) return -1;
  // Read and parse file contents
  fclose(fp);
  return 0;
}
\endcode

@subsection subsec_expert_step4 Step 4: Register Reader
Add call to ReadDefFileIdxPara() in readdef.c:
\code{c}
if (ReadMyNewFile(cFileNames[KWMyNewFile], &X->Def) != 0) {
  return -1;
}
\endcode

@section sec_expert_format File Format Guidelines
- First line: Number of entries (header)
- Comment lines start with @c #
- Use space or tab as delimiters
- Indices are 0-based

@page page_addmodpara Add new parameter into modpara

@section sec_modpara_overview Overview

ModPara file contains calculation parameters.
Parameters are read in readdef.c function ReadModPara().

@section sec_modpara_steps Steps to Add New Parameter

@subsection subsec_modpara_step1 Step 1: Add to DefineList
In struct.h, add new member to DefineList:
\code{c}
struct DefineList {
  // ... existing members ...
  int MyNewParam;    // Description of new parameter
};
\endcode

@subsection subsec_modpara_step2 Step 2: Add Keyword
In readdef.c, add keyword handling:
\code{c}
else if (strcmp(ctmp, "mynewparam") == 0) {
  fscanf(fp, "%d", &(Def->MyNewParam));
}
\endcode

@subsection subsec_modpara_step3 Step 3: Set Default Value
In setmem_def() (xsetmem.c), set default:
\code{c}
X->Def.MyNewParam = 0;  // Default value
\endcode

@subsection subsec_modpara_step4 Step 4: Validate
Add validation in check.c if needed:
\code{c}
if (X->Def.MyNewParam < 0) {
  fprintf(stdoutMPI, "Error: MyNewParam must be >= 0\n");
  return MPIFALSE;
}
\endcode

@page page_addcalcmod Add new calculation mode into calcmod

@section sec_calcmod_overview Overview

Calculation modes are selected by @c CalcMod parameter.
Each mode has a dedicated CalcBy*.c file.

@section sec_calcmod_steps Steps to Add New Calculation Mode

@subsection subsec_calcmod_step1 Step 1: Define Mode Constant
In DefCommon.h:
\code{c}
#define MyNewCalc 10  // New calculation mode ID
\endcode

@subsection subsec_calcmod_step2 Step 2: Create Implementation
Create CalcByMyNew.c and CalcByMyNew.h:
\code{c}
// CalcByMyNew.h
#pragma once
#include "Common.h"
int CalcByMyNew(struct EDMainCalStruct *X);

// CalcByMyNew.c
#include "CalcByMyNew.h"
int CalcByMyNew(struct EDMainCalStruct *X) {
  // Implementation
  return TRUE;
}
\endcode

@subsection subsec_calcmod_step3 Step 3: Add to Main Switch
In HPhiMain.c, add case to switch statement:
\code{c}
case MyNewCalc:
  if (CalcByMyNew(&X) != TRUE) {
    exitMPI(-3);
  }
  break;
\endcode

@subsection subsec_calcmod_step4 Step 4: Register in CMake
Add source file to src/CMakeLists.txt SOURCES list.

@page page_time Compute elapsed time for new functions

@section sec_time_overview Overview

HPhi uses hierarchical timer IDs defined in CalcTime.h/c.
Timer results are output to @c output/Time_*.dat.

@section sec_time_usage Timer Usage

@subsection subsec_time_basic Basic Usage
\code{c}
#include "CalcTime.h"

StartTimer(MyTimerID);
// ... code to measure ...
StopTimer(MyTimerID);
\endcode

@subsection subsec_time_ids Timer ID Convention
Timer IDs are hierarchical (3-digit numbers):
- @c 0 : Total execution time
- @c 1000 : sz() - Hilbert space construction
- @c 2000 : diagonalcalc() - Diagonal elements
- @c 3000 : TPQ calculation
- @c 4000 : Lanczos calculation
- @c 5000 : Full diagonalization
- @c 6000 : Spectrum calculation

Sub-timers use additional digits:
- @c 300 : Transfer terms (total)
- @c 310 : Local transfer
- @c 311-313 : Local transfer details
- @c 320 : InterAll terms
- @c 321-322 : InterAll details

@subsection subsec_time_new Adding New Timer

Choose an unused ID following the hierarchy:
\code{c}
#define TimerMyNew 350  // Under 300 category

StartTimer(TimerMyNew);
// ... new code ...
StopTimer(TimerMyNew);
\endcode

@section sec_time_output Timer Output

Timer results are written by OutputTimer():
- File: @c output/Time_*.dat
- Format: Timer ID, elapsed time, description

@page page_setmem Malloc vectors

@section sec_setmem_overview Overview

Memory allocation in HPhi uses dedicated functions in setmemory.c.
Main allocation routines are in xsetmem.c.

@section sec_setmem_functions Allocation Functions

@subsection subsec_setmem_1d 1D Arrays
\code{c}
int *arr_i = i_1d_allocate(N);           // int[N]
double *arr_d = d_1d_allocate(N);        // double[N]
double complex *arr_c = cd_1d_allocate(N); // complex[N]
unsigned long int *arr_l = lui_1d_allocate(N); // ulong[N]
\endcode

@subsection subsec_setmem_2d 2D Arrays
\code{c}
int **arr2d = i_2d_allocate(N, M);       // int[N][M]
double **arr2d = d_2d_allocate(N, M);    // double[N][M]
double complex **arr2d = cd_2d_allocate(N, M); // complex[N][M]
\endcode

2D arrays use contiguous memory allocation for cache efficiency.

@subsection subsec_setmem_3d 3D Arrays
\code{c}
int ***arr3d = i_3d_allocate(N, M, L);   // int[N][M][L]
double complex ***arr3d = cd_3d_allocate(N, M, L);
\endcode

@section sec_setmem_free Deallocation

Use corresponding free functions:
\code{c}
free_i_1d_allocate(arr_i);
free_cd_2d_allocate(arr2d);
\endcode

@section sec_setmem_routines Main Allocation Routines

@subsection subsec_setmem_head setmem_HEAD()
Allocates BindStruct and basic arrays.
Called at program start.

@subsection subsec_setmem_def setmem_def()
Allocates DefineList arrays after reading input sizes:
- Transfer arrays
- Interaction arrays
- Site information

@subsection subsec_setmem_large setmem_large()
Allocates large wavefunction vectors:
- @c v0, @c v1 : Main vectors (size: idim_max)
- @c v1buf : MPI buffer
- @c list_1, @c list_2_1, @c list_2_2 : Basis lists

This is the largest memory consumer.

@section sec_setmem_estimate Memory Estimation

Total memory scales as:
\f[
M \approx 3 \times D \times 16 \text{ bytes}
\f]
where \f$D\f$ is the Hilbert space dimension.
For Hubbard model: \f$D = \binom{N_s}{N_\uparrow}\binom{N_s}{N_\downarrow}\f$
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

  if (ValidateSymmetryRuntimeOptions(&(X.Bind)) != 0) {
    exitMPI(-1);
  }

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
  if (ValidateSymmetryHamiltonian(&(X.Bind)) != 0) {
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

    if (X.Bind.Def.iFlgSymmetryBasis == TRUE) {
      StartTimer(1100);
      if (BuildSymmetryBasis(&(X.Bind)) != 0) {
        StopTimer(1100);
        exitMPI(-1);
      }
      if (ActivateSymmetryBasisDimension(&(X.Bind)) != 0) {
        StopTimer(1100);
        exitMPI(-1);
      }
      if (ValidateSymmetrySectorOptions(&(X.Bind)) != 0) {
        StopTimer(1100);
        exitMPI(-1);
      }
      StopTimer(1100);
      StartTimer(1101);
      if (BuildSymmetryMatvecPlan(&(X.Bind)) != 0) {
        StopTimer(1101);
        exitMPI(-1);
      }
      StopTimer(1101);
    }
      
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
          StopTimer(5000);
          exitMPI(-3);
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
  FreeSymmetryBasis(X.Bind.Sym);
  X.Bind.Sym = NULL;
  FinalizeMPI();
  return 0;
}

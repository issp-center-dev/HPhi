/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "HΦ", "index.html", [
    [ "Coding rule", "page_codingrule.html", [
      [ "General Rules", "page_codingrule.html#sec_coding_general", null ],
      [ "Naming Conventions", "page_codingrule.html#sec_coding_naming", [
        [ "Function Names", "page_codingrule.html#subsec_functions", null ],
        [ "Variable Names", "page_codingrule.html#subsec_variables", null ],
        [ "Macro Names", "page_codingrule.html#subsec_macros", null ]
      ] ],
      [ "OpenMP Parallelization", "page_codingrule.html#sec_coding_openmp", null ],
      [ "MPI Parallelization", "page_codingrule.html#sec_coding_mpi", null ],
      [ "Memory Management", "page_codingrule.html#sec_coding_memory", null ],
      [ "Error Handling", "page_codingrule.html#sec_coding_error", null ],
      [ "Testing", "page_codingrule.html#sec_coding_testing", null ]
    ] ],
    [ "Add new source-file, executable, scripts (handle CMake)", "page_cmake.html", [
      [ "New source code", "page_cmake.html#sec_newsource", null ],
      [ "New executable", "page_cmake.html#sec_newexecutable", null ],
      [ "New script", "page_cmake.html#sec_newscript", null ],
      [ "CMake Build Options", "page_cmake.html#sec_cmake_options", null ]
    ] ],
    [ "Global variables and Data structure", "page_variable.html", [
      [ "Global Variables", "page_variable.html#sec_global_vars", [
        [ "MPI-related Variables", "page_variable.html#subsec_mpi_vars", null ],
        [ "Path and File Variables", "page_variable.html#subsec_path_vars", null ],
        [ "Wavefunction Vectors", "page_variable.html#subsec_wave_vars", null ]
      ] ],
      [ "Main Data Structures", "page_variable.html#sec_data_struct", [
        [ "BindStruct (X)", "page_variable.html#subsec_bindstruct", null ],
        [ "DefineList (X->Def)", "page_variable.html#subsec_definelist", null ],
        [ "LargeList (X->Large)", "page_variable.html#subsec_largelist", null ],
        [ "PhysList (X->Phys)", "page_variable.html#subsec_physlist", null ]
      ] ]
    ] ],
    [ "Various kind of Logs", "page_log.html", [
      [ "Progress Messages", "page_log.html#sec_log_progress", null ],
      [ "Error Messages", "page_log.html#sec_log_error", null ],
      [ "Time Logging", "page_log.html#sec_log_time", null ],
      [ "Debug Output", "page_log.html#sec_log_debug", null ]
    ] ],
    [ "Add new input-file for Expert mode", "page_addexpert.html", [
      [ "Overview", "page_addexpert.html#sec_expert_overview", null ],
      [ "Steps to Add New Input File", "page_addexpert.html#sec_expert_steps", [
        [ "Step 1: Define Keywords", "page_addexpert.html#subsec_expert_step1", null ],
        [ "Step 2: Add File Type", "page_addexpert.html#subsec_expert_step2", null ],
        [ "Step 3: Implement Reader", "page_addexpert.html#subsec_expert_step3", null ],
        [ "Step 4: Register Reader", "page_addexpert.html#subsec_expert_step4", null ]
      ] ],
      [ "File Format Guidelines", "page_addexpert.html#sec_expert_format", null ]
    ] ],
    [ "Add new parameter into modpara", "page_addmodpara.html", [
      [ "Overview", "page_addmodpara.html#sec_modpara_overview", null ],
      [ "Steps to Add New Parameter", "page_addmodpara.html#sec_modpara_steps", [
        [ "Step 1: Add to DefineList", "page_addmodpara.html#subsec_modpara_step1", null ],
        [ "Step 2: Add Keyword", "page_addmodpara.html#subsec_modpara_step2", null ],
        [ "Step 3: Set Default Value", "page_addmodpara.html#subsec_modpara_step3", null ],
        [ "Step 4: Validate", "page_addmodpara.html#subsec_modpara_step4", null ]
      ] ]
    ] ],
    [ "Add new calculation mode into calcmod", "page_addcalcmod.html", [
      [ "Overview", "page_addcalcmod.html#sec_calcmod_overview", null ],
      [ "Steps to Add New Calculation Mode", "page_addcalcmod.html#sec_calcmod_steps", [
        [ "Step 1: Define Mode Constant", "page_addcalcmod.html#subsec_calcmod_step1", null ],
        [ "Step 2: Create Implementation", "page_addcalcmod.html#subsec_calcmod_step2", null ],
        [ "Step 3: Add to Main Switch", "page_addcalcmod.html#subsec_calcmod_step3", null ],
        [ "Step 4: Register in CMake", "page_addcalcmod.html#subsec_calcmod_step4", null ]
      ] ]
    ] ],
    [ "Compute elapsed time for new functions", "page_time.html", [
      [ "Overview", "page_time.html#sec_time_overview", null ],
      [ "Timer Usage", "page_time.html#sec_time_usage", [
        [ "Basic Usage", "page_time.html#subsec_time_basic", null ],
        [ "Timer ID Convention", "page_time.html#subsec_time_ids", null ],
        [ "Adding New Timer", "page_time.html#subsec_time_new", null ]
      ] ],
      [ "Timer Output", "page_time.html#sec_time_output", null ]
    ] ],
    [ "Malloc vectors", "page_setmem.html", [
      [ "Overview", "page_setmem.html#sec_setmem_overview", null ],
      [ "Allocation Functions", "page_setmem.html#sec_setmem_functions", [
        [ "1D Arrays", "page_setmem.html#subsec_setmem_1d", null ],
        [ "2D Arrays", "page_setmem.html#subsec_setmem_2d", null ],
        [ "3D Arrays", "page_setmem.html#subsec_setmem_3d", null ]
      ] ],
      [ "Deallocation", "page_setmem.html#sec_setmem_free", null ],
      [ "Main Allocation Routines", "page_setmem.html#sec_setmem_routines", [
        [ "setmem_HEAD()", "page_setmem.html#subsec_setmem_head", null ],
        [ "setmem_def()", "page_setmem.html#subsec_setmem_def", null ],
        [ "setmem_large()", "page_setmem.html#subsec_setmem_large", null ]
      ] ],
      [ "Memory Estimation", "page_setmem.html#sec_setmem_estimate", null ]
    ] ],
    [ "Data Structures", "annotated.html", [
      [ "Data Structures", "annotated.html", "annotated_dup" ],
      [ "Data Structure Index", "classes.html", null ],
      [ "Data Fields", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Variables", "functions_vars.html", "functions_vars" ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "Globals", "globals.html", [
        [ "All", "globals.html", "globals_dup" ],
        [ "Functions", "globals_func.html", "globals_func" ],
        [ "Variables", "globals_vars.html", "globals_vars" ],
        [ "Typedefs", "globals_type.html", null ],
        [ "Enumerations", "globals_enum.html", null ],
        [ "Enumerator", "globals_eval.html", null ],
        [ "Macros", "globals_defs.html", "globals_defs" ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"CG__EigenVector_8c.html",
"ErrorMessage_8c.html#a555e61e3812f4a7501983309bf343d36",
"LogMessage_8c.html#a68e859f29dcdaa3a7a1c4f00b3885757",
"ProgressMessage_8c.html#aead2c0d63388c7b0fbdadf9345bbbd09",
"diagonalcalc_8c.html#a8aa25ee91e1e969763cdfdcb95d50456",
"functions_vars.html",
"global_8h.html#a9705c5a555c717a9667442bd8f6f2c6b",
"index.html",
"mltplyMPIBoost_8c.html",
"mltplySpinCore_8h.html#a51041879d0e2877713d258bdb75a4a3a",
"output_8h.html",
"restart__checkpoint_8h.html#a358fe4cfa60979d554eeae58643cb634",
"structCorrelationOrbit.html#a97cd870e15e221e6d53a432a8a704588",
"structGreenOutputManifestRecord.html#a1e64ffc6533b6e10d251fa26d7f6e2e1",
"structSymmetryBasisDistributionStats.html#a9633446683109bf786c5329583a6a261",
"structSymmetryMpiExchangeStats.html#ad1c148055fabeb03569b7c68ec0bde4c",
"symmetry__basis_8c.html#abc7c0ae431b0f2587cd4596948030d36",
"symmetry__distribution_8h.html#a5b824ae8652e0560ecd32683d1e4c02a",
"symmetry__mpi__exchange_8c.html#ad8a5c99900f5a0646121589fcbe7e3dc",
"version__major_8h_source.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';
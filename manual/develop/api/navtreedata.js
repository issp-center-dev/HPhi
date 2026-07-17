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
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"CG__EigenVector_8c.html",
"ErrorMessage_8c.html#ada7d332136387c8136cdfa9f35bd4536",
"LogMessage_8c.html#af71fc60a8638ffb4b8aa42613dbae0dc",
"anomalous__pair_8c.html#a554352de6619d8d7d94465a0f0ba6f4e",
"expec__cisajscktaltdc_8c.html#a5780ffd65b489119e7d7eb9617e2d3a0",
"global_8h.html#a18b86ba12763e11b6cbc878e43edb711",
"input_8c.html#a53da3ac7824590302b0aec6d7c9ecd57",
"mltplyMPIHubbardCore_8h.html#a70d6e1c992ec878c2c8ec3e8128e21ed",
"mltplySpin_8h.html#a367fcb99b257491d4370cdef55e496de",
"page_log.html#sec_log_debug",
"structDefineList.html#a0a1a6ed4ebd530f2986b6f3bf6ee3191",
"structMPIInterAllGroup.html#a0615b287419b5428b9b7c07d0f35f16a",
"sz_8h.html#ac329cb57c6b95deb1189bad49505e9ba"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';
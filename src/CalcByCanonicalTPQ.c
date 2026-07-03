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
#include "FirstMultiply.h"
#include "Multiply.h"
#include "expec_energy_flct.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "nbody_correlation.h"
#include "anomalous_pair.h"
#include "green_output.h"
#include "MakeIniVec.h"
#include "CalcByCanonicalTPQ.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#include <ctype.h>
#ifdef MPI
    #include <mpi.h>
#endif

/**
 * @file   CalcByCanonicalTPQ.c 
 * @version 3.4
 * @author Takahiro Misawa (BAQIS)
 * @brief  File for giving functions of the canonical TPQ (cTPQ) method
 *
 */

/** 
 * 
 * @brief A main function to calculate physical quqntities by the cTPQ method
 *
 * @param [in] NumAve  Number of samples
 * @param [in] ExpecInterval interval steps between the steps to calculate physical quantities
 * @param [in,out] X CalcStruct list for getting and giving calculation information
 * 
 * @author Takahiro Misawa (BAQIS)
 *
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 */
int CalcByCanonicalTPQ(
      const int NumAve,
      const int ExpecInterval,
      struct EDMainCalStruct *X
)
{
    char sdt[D_FileNameMax];
    char sdt_phys[D_FileNameMax];
    char sdt_norm[D_FileNameMax];
    char sdt_flct[D_FileNameMax];
    char file_name[D_FileNameMax];
    int rand_i, rand_max, iret;
    unsigned long int i_max;
    int step_i,step_iO=0;
    FILE *fp;
    double inv_temp,Ns,delta_tau;
    struct TimeKeepStruct tstruct;
    size_t byte_size;
    int green_output_initialized = 0;
    /*[s] for inverse temperatures*/
    double *read_invtemp=NULL;
    int    *read_nmax=NULL,*read_physcal=NULL,*read_eigen=NULL; 
    int    flag_read_invtemp;
    int    num_lines=0, read_lines=0, invtemp_status=0;
    /*[e] for inverse temperatures*/

    tstruct.tstart=time(NULL);
    
    rand_max  = NumAve;
    int step_spin = ExpecInterval;

    flag_read_invtemp = X->Bind.Def.flag_read_invtemp;
    /*[s]Following copilot's suggestion, we use strncpy to avoid buffer overflow*/
    strncpy(file_name, X->Bind.Def.file_invtemp, D_FileNameMax - 1);
    file_name[D_FileNameMax - 1] = '\0'; // Ensure null termination
    /*[e]Following copilot's suggestion, we use strncpy to avoid buffer overflow*/
    if (X->Bind.Def.flag_read_invtemp==1){
        if(X->Bind.Def.iReStart==RESTART_INOUT || X->Bind.Def.iReStart==RESTART_IN){
            fprintf(stdoutMPI,
                    "Error: InvTemp input cannot be combined with Restart=2 or Restart=3 in cTPQ.\n");
            return -1;
        }
        if(myrank==0){
            num_lines      = count_file_lines(file_name); /*count lines of files*/
            if(num_lines <= 0){
                fprintf(stdoutMPI,
                        "Error: InvTemp file '%s' must contain at least one complete row.\n",
                        file_name);
                invtemp_status = -1;
            }else{
                //[s] allocate
                read_invtemp   = (double*)calloc(num_lines,sizeof(double));
                read_nmax      = (int *)calloc(num_lines, sizeof(int));
                read_physcal   = (int *)calloc(num_lines, sizeof(int));
                read_eigen     = (int *)calloc(num_lines, sizeof(int));
                //[e] allocate
                if(read_invtemp == NULL || read_nmax == NULL ||
                   read_physcal == NULL || read_eigen == NULL){
                    fprintf(stdoutMPI,
                            "Error: failed to allocate InvTemp arrays for %d rows.\n",
                            num_lines);
                    invtemp_status = -1;
                }else{
                    read_lines = func_read_invtemp(read_invtemp,read_nmax,read_physcal,read_eigen,file_name,num_lines); /*read files*/
                    if(read_lines < 0){
                        fprintf(stdoutMPI,
                                "Error: InvTemp file '%s' has an invalid row. Each non-empty row must contain exactly: beta nmax physcal eigen.\n",
                                file_name);
                        invtemp_status = -1;
                    }else if(read_lines == 0){
                        fprintf(stdoutMPI,
                                "Error: InvTemp file '%s' must contain at least one complete row.\n",
                                file_name);
                        invtemp_status = -1;
                    }else{
                        num_lines = read_lines;
                    }
                }
            }
        }
        #ifdef MPI
            MPI_Bcast(&invtemp_status, 1, MPI_INT, 0, MPI_COMM_WORLD);
            // Broadcast the number of lines to all ranks
            MPI_Bcast(&num_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
            //printf("DEBUG: myrank = %d: flag_read_invtemp = %d file_name = %s num_lines = %d \n", 
            //       myrank,X->Bind.Def.flag_read_invtemp, file_name, num_lines);
        #endif

        if(invtemp_status != 0){
            free(read_invtemp);
            free(read_nmax);
            free(read_physcal);
            free(read_eigen);
            return -1;
        }

        #ifdef MPI
            int local_invtemp_status = 0;
            // Allocate memory on non-root ranks
            if (myrank != 0) {
                read_invtemp = (double*)calloc(num_lines, sizeof(double));
                read_nmax    = (int*)calloc(num_lines, sizeof(int));
                read_physcal = (int*)calloc(num_lines, sizeof(int));
                read_eigen   = (int*)calloc(num_lines, sizeof(int));
                if(read_invtemp == NULL || read_nmax == NULL ||
                   read_physcal == NULL || read_eigen == NULL){
                    fprintf(stdoutMPI,
                            "Error: failed to allocate InvTemp arrays for %d rows.\n",
                            num_lines);
                    local_invtemp_status = -1;
                }
            }
            MPI_Allreduce(&local_invtemp_status, &invtemp_status, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            if(invtemp_status != 0){
                free(read_invtemp);
                free(read_nmax);
                free(read_physcal);
                free(read_eigen);
                return -1;
            }

            // Broadcast data arrays to all ranks
            MPI_Bcast(read_invtemp, num_lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(read_nmax,    num_lines, MPI_INT,    0, MPI_COMM_WORLD);
            MPI_Bcast(read_physcal, num_lines, MPI_INT,    0, MPI_COMM_WORLD);
            MPI_Bcast(read_eigen,   num_lines, MPI_INT,    0, MPI_COMM_WORLD);        
        #endif
    }else{
        if (X->Bind.Def.Param.ExpandCoef==0){
            X->Bind.Def.Param.ExpandCoef=10;   
            fprintf(stdoutMPI, "In cTPQ calc., the default value of ExpandCoef (=10) is used. \n");
        }else{
            fprintf(stdoutMPI, "In cTPQ calc., ExpandCoef is specified as %d. \n",X->Bind.Def.Param.ExpandCoef);
        }
    }
    X->Bind.Def.St=0;
    fprintf(stdoutMPI, "%s", cLogTPQ_Start);
    for (rand_i = 0; rand_i<rand_max; rand_i++){
        if(X->Bind.Def.iOutputDataHead==1){
            int prefix_length;
            prefix_length = sprintf(sdt_phys, "%s_", X->Bind.Def.CDataFileHead);
            sprintf(sdt_phys + prefix_length, cFileNameSSRand, rand_i);
            prefix_length = sprintf(sdt_norm, "%s_", X->Bind.Def.CDataFileHead);
            sprintf(sdt_norm + prefix_length, cFileNameNormRand, rand_i);
            prefix_length = sprintf(sdt_flct, "%s_", X->Bind.Def.CDataFileHead);
            sprintf(sdt_flct + prefix_length, cFileNameFlctRand, rand_i);
        }else{
            sprintf(sdt_phys, cFileNameSSRand, rand_i);
            sprintf(sdt_norm, cFileNameNormRand, rand_i);
            sprintf(sdt_flct, cFileNameFlctRand, rand_i);
        }
        Ns = 1.0 * X->Bind.Def.NsiteMPI;
        fprintf(stdoutMPI, cLogTPQRand, rand_i+1, rand_max);
        iret=0;
        X->Bind.Def.irand=rand_i;

        //Make or Read initial vector
        if(X->Bind.Def.iReStart==RESTART_INOUT || X->Bind.Def.iReStart==RESTART_IN) {
            StartTimer(3600);
            TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
            fprintf(stdoutMPI, "%s", cLogInputVecStart);
            sprintf(sdt, cFileNameInputVector, rand_i, myrank);
            childfopenALL(sdt, "rb", &fp);
            if(fp==NULL){
                fprintf(stderr, "A file of Inputvector does not exist (rank %d).\n", myrank);
                iret=1;
            }
            /* All ranks must agree on the fallback: tmpvec files are per-rank and a
               partially missing set would otherwise split ranks across the restart
               and fresh-start branches, whose MPI collectives do not match. */
            iret = (int)MaxMPI_li((unsigned long int)iret);
            if(iret==1){
                if(fp != NULL) fclose(fp);
                fprintf(stdoutMPI, "Start to calculate in normal procedure.\n");
                StopTimer(3600);
            }else{
                byte_size = fread(&step_i, sizeof(step_i), 1, fp);
                byte_size = fread(&i_max, sizeof(long int), 1, fp);
                if(i_max != X->Bind.Check.idim_max){
                    fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
                    exitMPI(-1);
                }
                byte_size = fread(v0, sizeof(complex double), X->Bind.Check.idim_max+1, fp);
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
                fprintf(stdoutMPI, "%s", cLogInputVecFinish);
                fclose(fp);
                StopTimer(3600);
                X->Bind.Def.istep=step_i;
                StartTimer(3200);
                iret=expec_energy_flct(&(X->Bind)); //v1 <- v0 and v0 = H*v1
                StopTimer(3200);
                if(iret != 0) return -1;
                step_iO = step_i-1;
                if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
            }
        }
        if(X->Bind.Def.iReStart==RESTART_NOT || X->Bind.Def.iReStart==RESTART_OUT || iret ==1) {
            StartTimer(3600);
            if (green_output_initialized == 0) {
                if (GreenOutputInitializeAggregateFiles(&(X->Bind)) != 0) {
                    return -1;
                }
                green_output_initialized = 1;
            }
            if (childfopenMPI(sdt_phys, "w", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%s", cLogSSRand);
            fclose(fp);
            // for norm
            if (childfopenMPI(sdt_norm, "w", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%s", cLogNormRand);
            fclose(fp);
            // for fluctuations
            if (childfopenMPI(sdt_flct, "w", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%s", cLogFlctRand);
            fclose(fp);
            StopTimer(3600);
            step_i = 0;
            StartTimer(3100);
            if(rand_i==0){
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "w", rand_i, step_i);
            }else{
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
            }
            /**@brief
             Initialize v1 and v0 = v1
            */
            MakeIniVec(rand_i, &(X->Bind)); 
            /*[s] tau*/
            inv_temp  = 0.0;
            delta_tau = 1.0/LargeValue;
            /*[e] tau*/
            StopTimer(3100);
            // for norm
            if (childfopenMPI(sdt_norm, "a", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, global_1st_norm, global_1st_norm, step_i);
            fclose(fp);
            /**@brief
             Compute expectation value at infinite temperature
            */
            X->Bind.Def.istep = 0;
            StartTimer(3300);
            iret=expec_cisajs(&(X->Bind), v1);
            StopTimer(3300);
            if(iret !=0) return -1;

            StartTimer(3400);
            iret=expec_cisajscktaltdc(&(X->Bind), v1);
            StopTimer(3400);
            if(iret !=0) return -1;
            iret=expec_nbodyg(&(X->Bind), v1);
            if(iret !=0) return -1;
            iret=expec_anomalousg(&(X->Bind), v1);
            if(iret !=0) return -1;
        
            StartTimer(3200);
            iret=expec_energy_flct(&(X->Bind)); //v1 <- v0 and v0 = H*v1
            StopTimer(3200);
            if(iret !=0) return -1;
            //inv_temp = 0; /* (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);*/
            if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var,
            X->Bind.Phys.doublon, X->Bind.Phys.num, step_i);
            fclose(fp);
            StartTimer(3600);
            // for fluctuations
            if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp,X->Bind.Phys.num,X->Bind.Phys.num2, X->Bind.Phys.doublon,X->Bind.Phys.doublon2, X->Bind.Phys.Sz,X->Bind.Phys.Sz2,step_i);
            fclose(fp);
            StopTimer(3600);
            step_i += 1;
            X->Bind.Def.istep = step_i;
            step_iO=0;
        }
        if (flag_read_invtemp == 1){
            X->Bind.Def.Lanczos_max = num_lines;
            if(read_eigen[0]==1){
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, 0);
                fprintf(stdoutMPI, "%s", cLogOutputVecStart);
                sprintf(sdt, cFileNameTPQVector, rand_i, myrank,0);
                if(childfopenALL(sdt, "wb", &fp)!=0){
                    exitMPI(-1);
                }
                fwrite(&step_i, sizeof(step_i), 1, fp);
                fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
                fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
                fclose(fp);
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, 0);
                fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
            } 
            //printf("num_lines = %d  %d %d \n",num_lines,X->Bind.Def.istep,X->Bind.Def.Lanczos_max);
        }
        int progress_interval = (X->Bind.Def.Lanczos_max - step_iO) / 10;
        if (progress_interval < 1) progress_interval = 1;  /* avoid % 0 for small (Lanczos_max - step_iO) */
        for (step_i = X->Bind.Def.istep; step_i<X->Bind.Def.Lanczos_max; step_i++){
            if (flag_read_invtemp == 1){
                X->Bind.Def.Param.ExpandCoef = read_nmax[step_i-1];
                delta_tau                    = read_invtemp[step_i]-read_invtemp[step_i-1];
            }
            if(step_i % progress_interval == 0){
                fprintf(stdoutMPI, cLogTPQStep, step_i, X->Bind.Def.Lanczos_max);
            }
            X->Bind.Def.istep=step_i;
            StartTimer(3600);
            TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
            StopTimer(3600);
            StartTimer(3500);
            MultiplyForCanonicalTPQ(&(X->Bind),delta_tau); // v0=exp[-delta_tau*H/2]*v1 in 4th order
            StopTimer(3500);

            StartTimer(3200);
            iret=expec_energy_flct(&(X->Bind)); //v1 <- v0 and v0 = H*v1
            StopTimer(3200);
            if(iret !=0) return -1;
            //inv_temp = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
            //printf("step_i = %d  inv_temp = %lf delta_tau = %lf\n", step_i, inv_temp,delta_tau);
            inv_temp  += delta_tau;
            //temp      = 1.0/inv_temp;

            StartTimer(3600);
            if(childfopenMPI(sdt_phys, "a", &fp)!=0){
                return FALSE;
            }
            fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var, X->Bind.Phys.doublon, X->Bind.Phys.num ,step_i);
            fclose(fp);

            if(childfopenMPI(sdt_norm, "a", &fp)!=0){
                return FALSE;
            }
            fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, global_norm, global_1st_norm, step_i);
            fclose(fp);

            // for fluctuations
            if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
                return -1;
            }
            fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp,X->Bind.Phys.num,X->Bind.Phys.num2, X->Bind.Phys.doublon,X->Bind.Phys.doublon2, X->Bind.Phys.Sz,X->Bind.Phys.Sz2,step_i);
            fclose(fp);
            StopTimer(3600);

            if (flag_read_invtemp == 1){
                if (read_physcal[step_i] == 1){
                    StartTimer(3300);
                    iret=expec_cisajs(&(X->Bind),v1);
                    StopTimer(3300);
                    if(iret !=0) return -1;

                    StartTimer(3400);
                    iret=expec_cisajscktaltdc(&(X->Bind), v1);
                    StopTimer(3400);
                    if(iret !=0) return -1;
                    iret=expec_nbodyg(&(X->Bind), v1);
                    if(iret !=0) return -1;
                    iret=expec_anomalousg(&(X->Bind), v1);
                    if(iret !=0) return -1;
                }
            }else{  
                if (step_i%step_spin == 0){
                    StartTimer(3300);
                    iret=expec_cisajs(&(X->Bind),v1);
                    StopTimer(3300);
                    if(iret !=0) return -1;

                    StartTimer(3400);
                    iret=expec_cisajscktaltdc(&(X->Bind), v1);
                    StopTimer(3400);
                    if(iret !=0) return -1;
                    iret=expec_nbodyg(&(X->Bind), v1);
                    if(iret !=0) return -1;
                    iret=expec_anomalousg(&(X->Bind), v1);
                    if(iret !=0) return -1;
                }
            }
            if (flag_read_invtemp == 1 && read_eigen[step_i]==1){
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
                fprintf(stdoutMPI, "%s", cLogOutputVecStart);
                sprintf(sdt, cFileNameTPQVector, rand_i, myrank,step_i);
                if(childfopenALL(sdt, "wb", &fp)!=0){
                    exitMPI(-1);
                }
                fwrite(&step_i, sizeof(step_i), 1, fp);
                fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
                fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
                fclose(fp);
                TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
                fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
            }  
        }
        if(X->Bind.Def.iReStart== RESTART_OUT || X->Bind.Def.iReStart==RESTART_INOUT){
            TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
            fprintf(stdoutMPI, "%s", cLogOutputVecStart);
            sprintf(sdt, cFileNameOutputVector, rand_i, myrank);
            if(childfopenALL(sdt, "wb", &fp)!=0){
                exitMPI(-1);
            }
            fwrite(&step_i, sizeof(step_i), 1, fp);
            fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
            fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
            fclose(fp);
            TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
            fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
        }
    }
    fprintf(stdoutMPI, "%s", cLogTPQ_End);

    tstruct.tend=time(NULL);
    fprintf(stdoutMPI, cLogTPQEnd, (int)(tstruct.tend-tstruct.tstart));
    /*[s] Free memory for inverse temperature data if it was read from a file */
    if (flag_read_invtemp == 1){
        free(read_invtemp);
        free(read_nmax);
        free(read_physcal);
        free(read_eigen);
    }
    /*[e] Free memory for inverse temperature data if it was read from a file */
    return TRUE;
}

/**
 * @brief Count the number of lines in a file
 *
 * @param file_name Name of the file to count lines in
 * @return Number of lines in the file, or -1 if the file could not be opened
 * @author Takahiro Misawa (ISSP)
 */
int count_file_lines(const char *file_name) {
    FILE *file;
    int lines = 0;
    int ch;
    int prev = '\n'; // assume file starts with a new line

    file = fopen(file_name, "r");
    if (file == NULL) {
        fprintf(stderr, "could not open file: %s\n", file_name);
        return -1; // distinguishable error code
    }

    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            lines++;
        }
        prev = ch;
    }

    fclose(file);

    if (prev != '\n') {
        lines++;
    }

    return lines;
}

/* * @brief Read inverse temperature data from a file
 *
 * @param read_invtemp Pointer to an array for storing inverse temperatures
 * @param read_nmax Pointer to an array for storing nmax values
 * @param read_physcal Pointer to an array for storing physical calculation flags
 * @param read_eigen Pointer to an array for storing eigenvalue flags
 * @param file_name Name of the file to read from
 * @param max_lines Maximum number of lines to read from the file
 * @return Number of lines read, or -1 if the file could not be opened, or -2 if the number of lines exceeds max_lines
 * @author Takahiro Misawa (ISSP)
 */
int func_read_invtemp(double *read_invtemp, int *read_nmax, int *read_physcal, int *read_eigen, const char *file_name, int max_lines) {

    FILE *file = fopen(file_name, "r");
    char line[1024];
    if (file == NULL) {
        fprintf(stderr, "could not open file: %s\n", file_name);
        return -1;
    }

    int i = 0;
    int line_no = 0;
    while (fgets(line, sizeof(line), file) != NULL) {
        double invtemp_tmp;
        int nmax_tmp, physcal_tmp, eigen_tmp;
        int nread = 0;
        char *cursor = line;
        line_no++;

        while (isspace((unsigned char)*cursor)) cursor++;
        if (*cursor == '\0') continue;

        if (i >= max_lines) {
            fprintf(stderr, "Error: number of lines in file exceeds expected %d\n", max_lines);
            fclose(file);
            return -2;
        }

        if (sscanf(cursor, "%lf %d %d %d %n",
                   &invtemp_tmp, &nmax_tmp, &physcal_tmp, &eigen_tmp, &nread) != 4) {
            fprintf(stderr, "Error: invalid InvTemp row %d in file: %s\n", line_no, file_name);
            fclose(file);
            return -3;
        }
        cursor += nread;
        while (isspace((unsigned char)*cursor)) cursor++;
        if (*cursor != '\0') {
            fprintf(stderr, "Error: invalid InvTemp row %d in file: %s\n", line_no, file_name);
            fclose(file);
            return -3;
        }

        read_invtemp[i] = invtemp_tmp;
        read_nmax[i] = nmax_tmp;
        read_physcal[i] = physcal_tmp;
        read_eigen[i] = eigen_tmp;
        i++;
    }

    fclose(file);
    return i;
}

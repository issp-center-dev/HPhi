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
#include "matrixscalapack.h"
/**
 * @file matrixscalapack.c
 * @version 3.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 * 
 * @brief File for diagonalization using scalapack with functions for getting indices of global and local array
 * 
 *
 */

#ifdef _SCALAPACK
int use_scalapack = 0;

/**
 * @brief compute block size for scalapack
 * @param[in] Msize size of matrix (Msize x Msize)
 * @param[in] nproc number of processes
 * @return block size for scalapack
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetBlockSize(long int Msize, long int nproc) {
  long int block_size = 16;
  if(Msize*Msize/nproc > block_size*block_size)
    return block_size;
  return 1;
}

/**
 * @brief get processor array index from global array index
 * @param[in] i global array index
 * @param[in] np processor array dimension
 * @param[in] nb block size
 * @return processor array index
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetPArrayIndex(long int i, long int np, long int nb) {
  return (i/nb)%np;
}

/**
 * @brief get local array index from global array index
 * @param[in] i global array index
 * @param[in] np processor array dimension
 * @param[in] nb block size
 * @return local array index
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetLocalIndex(long int i, long int np, long int nb) {
  return (i/(np*nb))*nb + i%nb;
}

/**
 * @brief get global array index from local array index and processor array index
 * @param[in] il local array index
 * @param[in] p processor array index
 * @param[in] np processor array dimension
 * @param[in] nb block size
 * @return global array index
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetGlobalIndex(long int il, long int p, long int np, long int nb){
  return ((il/nb)*np+p)*nb + il%nb;
}

/**
 * @brief get rank of processor from indices of global matrix
 * @param[in] i index of global matrix
 * @param[in] j index of global matrix
 * @param[in] nprow processor array dimension for row
 * @param[in] npcol processor array dimension for column
 * @param[in] nb block size
 * @return rank of processor
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int MatToRank(long int i, long int j, long int nprow, long int npcol, long int nb){
  long int iproc, jproc;
  iproc = GetPArrayIndex(i, nprow, nb);
  jproc = GetPArrayIndex(j, npcol, nb);
  return iproc+jproc*nprow;
}

/**
 * @brief get column index of global matrix from given rank of processor and column index of local matrix
 * @param[in] lj column index of local matrix
 * @param[in] rank rank of processor
 * @param[in] npcol processor array dimension for column
 * @param[in] nb block size
 * @return column index of global matrix
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int GetMatRawInRank(long int lj, long int rank, long int npcol, long int nb){
  long int pcol = rank/npcol;
  return GetGlobalIndex(lj, pcol, npcol, nb);
}

/**
 * @brief get indices of local matrix from indices (need to free memory after this function used)
 * @param[in] i index of local matrix
 * @param[in] j index of local matrix
 * @param[in] nprow processor array dimension for row
 * @param[in] npcol processor array dimension for column
 * @param[in] nb block size
 * @return indices array of local matrix
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
long int *GetMatElementInRank(long int i, long int j, long int nprow, long int npcol, long int nb){
  long int *ij;
  ij = malloc(2*sizeof(int));
  ij[0] = GetLocalIndex(i, nprow, nb);
  ij[1] = GetLocalIndex(j, npcol, nb);
  return ij;
}

/**
 * @brief divide matrix
 * @param[in] m index of column of matrix
 * @param[in] n index of row of matrix
 * @param[in] Aorgmn value of matrix in (i, j)
 * @param[in, out] A distribution matrix
 * @param[in] desca descriptor of distribution matrix A
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
void DivMat(long int m, long int n, double complex Aorgmn, double complex *A, int *desca){
  long int mp = m+1, np = n+1;
  pzelset_(A, &mp, &np, desca, &Aorgmn);
}

/**
 * @brief get eigenvector from distributed matrix
 * @param[in] i index of eigenvector
 * @param[in] m size of eigenvector
 * @param[in] Z distribution matrix of eigenvector
 * @param[in] descZ descriptor for Z
 * @param[in, out] vec eigenvector
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
void GetEigenVector(long int i, long int m, double complex *Z, int *descZ, double complex *vec) {
  double complex alpha;
  long int j, ip, jp;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ip = i+1;
  for(j=0; j<m; j++){
    jp = j+1;
    pzelget_("A", " ", &alpha, Z, &jp, &ip, descZ);
    if(rank==0) {
      vec[j] = alpha;
    }
  }
}

/**
 * @brief diagonalization using scalapack
 * @param[in] xNsize size of matrix
 * @param[in] A input matrix
 * @param[in, out] r eigenvalue
 * @param[in, out] Z distribution matrix of eigenvector
 * @param[in, out] descZ descriptor for Z
 * @return this returns 0 when it finished normally
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
int diag_scalapack_cmp(long int xNsize, double complex **A, 
                       double complex *r, double complex *Z, int *descZ) {
  const int i_one=1, i_zero=0;
  const long int i_negone=-1;
  const double zero=0.0, one=1.0;
  long int m, n, mb, nb;
  int nprow, npcol;
  int myrow, mycol, info, lld;
  long int mp, nq;
  int ictxt;
  complex double *A_distr, *work, *rwork;
  double *W;
  int descA_distr[9];
  int rank, size, iam, nprocs;
  long int lwork, lrwork;
  int dims[2]={0,0};
  long int i, j, ip, jp;
  m=n=xNsize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size,2,dims);
  nprow=dims[0]; npcol=dims[1];
 
  blacs_pinfo_(&iam, &nprocs); 
  blacs_get_((int *)&i_negone, &i_zero, &ictxt);
  blacs_gridinit_(&ictxt, "R", &nprow, &npcol);
  blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);
 
  mb = GetBlockSize(m, size);
  nb = GetBlockSize(n, size);

  mp = numroc_(&m, &mb, &myrow, &i_zero, &nprow);
  nq = numroc_(&n, &nb, &mycol, &i_zero, &npcol);
  W = malloc(n*sizeof(double));
  A_distr = malloc(mp*nq*sizeof(complex double));

  lld = (mp>0) ? mp : 1;
  descinit_(descA_distr, &m, &n, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
  descinit_(descZ, &m, &n, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      DivMat(i, j, A[i][j], A_distr, descA_distr);
    }
  }

  double complex wkopt, rwkopt;
  pzheev_("V", "U", &n, A_distr, &i_one, &i_one, descA_distr, W, Z, &i_one, &i_one, descZ, &wkopt, &i_negone, &rwkopt, &i_negone, &info);

  lwork = (long int)wkopt;
  lrwork = (long int)rwkopt;
  work = malloc(lwork*sizeof(complex double));
  rwork = malloc(lrwork*sizeof(complex double));


  pzheev_("V", "U", &n, A_distr, &i_one, &i_one, descA_distr, W, Z, &i_one, &i_one, descZ, work, &lwork, rwork, &lrwork, &info);

  if(rank == 0){
    for(i=0; i<n; i++){
      r[i] = W[i];
    }
  }

  free(A_distr);
  free(work);
  free(rwork);
  free(W);

  use_scalapack = 1;

  return 0;
}

/* Cached destination grid for GetEigenVectorBlock: a 1x1 BLACS grid
   containing rank 0 only, built once per run (design doc section 3). */
static int ictxt_gather = -100;   /* -100: not initialized */
static int desc_gather[9];

/**
 * @brief Initialize the rank-0-only destination grid and descriptor.
 * All ranks must call this (blacs_gridmap_ is collective). Follows the
 * p?gemr2d contract: ranks outside the destination grid keep
 * desc[CTXT_] = -1 while all other descriptor fields stay valid.
 */
static void InitEigenVectorGatherContext(long int xNsize) {
  int i_negone = -1, i_zero = 0;
  int imap[1] = {0};
  int ld = 1, np_gather = 1;
  int myrow_g, mycol_g, nprow_g, npcol_g;

  blacs_get_(&i_negone, &i_zero, &ictxt_gather);
  blacs_gridmap_(&ictxt_gather, imap, &ld, &np_gather, &np_gather);

  /* Fully initialize the descriptor on ALL ranks (some implementations
     inspect fields other than CTXT_), then mark non-participants. */
  desc_gather[0] = 1;               /* DTYPE_: dense */
  desc_gather[1] = ictxt_gather;    /* CTXT_ */
  desc_gather[2] = (int)xNsize;     /* M_ */
  desc_gather[3] = 1;               /* N_ */
  desc_gather[4] = (int)xNsize;     /* MB_ */
  desc_gather[5] = 1;               /* NB_ */
  desc_gather[6] = 0;               /* RSRC_ */
  desc_gather[7] = 0;               /* CSRC_ */
  desc_gather[8] = (int)xNsize;     /* LLD_ */

  blacs_gridinfo_(&ictxt_gather, &nprow_g, &npcol_g, &myrow_g, &mycol_g);
  if (myrow_g < 0) {
    desc_gather[1] = -1;            /* not in the destination grid */
  }
}

/**
 * @brief get eigenvector from distributed matrix via a single block
 * transfer (pzgemr2d_) instead of the element-wise pzelget_ loop used
 * by GetEigenVector.
 * @param[in] idx 0-based index of eigenvector (column of Z)
 * @param[in] xNsize size of eigenvector
 * @param[in] Z distribution matrix of eigenvector
 * @param[in] descZ descriptor for Z
 * @param[in, out] vec eigenvector; on rank 0, vec[0..xNsize-1] holds the
 * gathered eigenvector. On other ranks vec is a work buffer with
 * unspecified contents and must be non-NULL.
 * @return 0 on success.
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Yusuke Konishi (Academeia Co., Ltd.)
 */
int GetEigenVectorBlock(long int idx, long int xNsize,
                        double complex *Z, int *descZ,
                        double complex *vec) {
  const long int i_one = 1;
  long int icol = idx + 1;
  long int m = xNsize, n = 1;

  if (ictxt_gather == -100) {
    InitEigenVectorGatherContext(xNsize);
  }
  /* Last argument: a context containing the union of both grids
     = the all-rank 2D context of Z. */
  pzgemr2d_(&m, &n, Z, (long int *)&i_one, &icol, descZ,
            vec, (long int *)&i_one, (long int *)&i_one, desc_gather,
            &descZ[1]);
  return 0;
}

/**
 * @brief Release the cached destination grid (call after the phys loop).
 */
void FreeEigenVectorGatherContext(void) {
  int myrow_g, mycol_g, nprow_g, npcol_g;
  if (ictxt_gather == -100) return;
  blacs_gridinfo_(&ictxt_gather, &nprow_g, &npcol_g, &myrow_g, &mycol_g);
  if (myrow_g >= 0) {
    blacs_gridexit_(&ictxt_gather);
  }
  ictxt_gather = -100;
}

/**
 * @brief Release a distributed eigenvector matrix and its BLACS contexts.
 *
 * The source 2D grid stored in descZ must remain alive while phys() gathers
 * or redistributes eigenvectors. Call this only after the observable phase.
 * Pointer arguments make the helper usable by focused unit tests without
 * depending on HPhi's process-global storage.
 */
void FreeDistributedEigenvectors(double complex **Z, int *descZ, int *used) {
  int ictxt;
  if (Z == NULL || descZ == NULL || used == NULL || !*used) return;
  FreeEigenVectorGatherContext();
  free(*Z);
  *Z = NULL;
  ictxt = descZ[1];
  if (ictxt >= 0) blacs_gridexit_(&ictxt);
  descZ[1] = -1;
  *used = 0;
}

/**
 * @brief Redistribute a 1D column-block panel into an existing 2D
 * block-cyclic matrix with a single pzgemr2d_ call (design doc sec. 3
 * phase 2). The 1D source: a 1 x P grid ('R'), MB = N (all rows in one
 * block), NB = NC = ceil(N/P) (one column block per rank), RSRC=CSRC=0,
 * LLD = panel_ld. Ranks owning zero columns still participate (their
 * local part is an unused 1-element buffer); pzgemr2d_ only touches the
 * columns numroc() assigns to that rank, so the buffer is never read.
 * @param[in] xNsize global matrix dimension N
 * @param[in] jbegin first owned column, 1-based (unused directly here;
 * kept for interface symmetry with the panel globals)
 * @param[in] ncols_panel number of columns owned by this rank (unused
 * directly here; ownership is re-derived from the 1D grid/NB so it must
 * match jbegin/ncols_panel by construction -- see xsetmem.c)
 * @param[in] panel_ld leading dimension of panel (= xNsize)
 * @param[in, out] panel this rank's 1D column-panel buffer
 * @param[in, out] A_distr destination 2D block-cyclic matrix
 * @param[in] descA_2d descriptor of A_distr
 * @return 0 on success, -1 on failure (same value on all ranks)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int RedistPanelToBlockCyclic(long int xNsize, long int jbegin,
                             long int ncols_panel, long int panel_ld,
                             double complex *panel,
                             double complex *A_distr, int *descA_2d) {
  int i_negone = -1, i_zero_i = 0, info;
  int ictxt_1d, nprow_1, npcol_1, myrow_1, mycol_1;
  int desc1d[9];
  int lld;
  long int NC, mb1, nb1;
  const long int i_one = 1;
  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  NC = (xNsize + size - 1) / size;

  blacs_get_(&i_negone, &i_zero_i, &ictxt_1d);
  nprow_1 = 1; npcol_1 = size;
  blacs_gridinit_(&ictxt_1d, "R", &nprow_1, &npcol_1);
  blacs_gridinfo_(&ictxt_1d, &nprow_1, &npcol_1, &myrow_1, &mycol_1);

  /* The 1x P 'R' grid must map mycol to this rank 1:1, matching
     xsetmem's panel ownership (jb = myrank*NC + 1). A mismatch would
     silently scramble columns, so check at runtime (survives NDEBUG)
     and synchronize the verdict so no rank enters the collective
     pzgemr2d_ alone (same pattern as SyncError in matrixlapack_elpa.c).
     Also verify that the caller-supplied ownership parameters
     (jbegin, ncols_panel) agree with the ownership derived internally
     from mycol_1/NC/xNsize -- this is what jbegin/ncols_panel are for,
     rather than being unused (void)-discarded. */
  {
    int ok = (mycol_1 == myrank) ? 0 : -1, gok;
    long int nc_expect = (mycol_1 < (int)((xNsize + NC - 1) / NC))
                           ? (((long int)mycol_1 + 1) * NC <= xNsize
                                ? NC : xNsize - (long int)mycol_1 * NC)
                           : 0;
    long int jb_expect = (long int)mycol_1 * NC + 1;
    if (ok != 0) {
      fprintf(stdout,
              "  Error: BLACS 1D grid column (%d) does not match MPI rank (%d):\n"
              "         panel ownership is inconsistent; aborting redistribution.\n",
              mycol_1, myrank);
    }
    if (ncols_panel != nc_expect ||
        (ncols_panel > 0 && jbegin != jb_expect)) {
      ok = -1;
      fprintf(stdout,
              "  Error: caller-supplied panel ownership (jbegin=%ld, ncols=%ld) does\n"
              "         not match the internally derived ownership (jbegin=%ld, ncols=%ld)\n"
              "         for rank %d: panel ownership is inconsistent; aborting redistribution.\n",
              jbegin, ncols_panel, jb_expect, nc_expect, myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      blacs_gridexit_(&ictxt_1d);
      return -1;
    }
  }

  lld = (panel_ld > 0) ? (int)panel_ld : 1;
  mb1 = xNsize;
  nb1 = NC;
  descinit_(desc1d, &xNsize, &xNsize, &mb1, &nb1, &i_zero_i, &i_zero_i,
            &ictxt_1d, &lld, &info);
  /* Synchronize the descinit_ verdict so no rank enters the collective
     pzgemr2d_ with an invalid descriptor (same pattern as
     RedistBlockCyclicToStatePanel below). */
  {
    int ok = (info == 0) ? 0 : -1, gok;
    if (ok != 0) {
      fprintf(stdout,
              "  Error: descinit_ failed (info=%d) for the generation-panel descriptor\n"
              "         on rank %d; aborting redistribution.\n", info, myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      blacs_gridexit_(&ictxt_1d);
      return -1;
    }
  }

  pzgemr2d_(&xNsize, &xNsize,
           panel, (long int *)&i_one, (long int *)&i_one, desc1d,
           A_distr, (long int *)&i_one, (long int *)&i_one, descA_2d,
           &descA_2d[1]);

  blacs_gridexit_(&ictxt_1d);
  return 0;
}

/**
 * @brief Redistribute an existing 2D block-cyclic matrix Z (e.g. the
 * eigenvector matrix produced by diag_scalapack_cmp/diag_elpa_cmp) into
 * a 1D state-column panel with a single pzgemr2d_ call -- the mirror
 * image of RedistPanelToBlockCyclic (source and destination roles
 * swapped; Z/descZ is the source here, panel is the destination). Same
 * 1 x P 'R' 1D grid, ownership formula, and zero-owner-participates
 * pattern as RedistPanelToBlockCyclic: NC = ceil(N/P),
 * first_state(r) = r*NC+1 (1-based), ncols_local(r) =
 * max(0, min((r+1)*NC, N) - r*NC). Ranks owning zero states still
 * participate with a 1-element dummy buffer + valid descriptor;
 * pzgemr2d_ only touches the columns numroc() assigns to that rank, so
 * the buffer is never read/written for those ranks.
 * @param[in] xNsize global matrix dimension N
 * @param[in] Z source 2D block-cyclic matrix (e.g. ScaLAPACK/ELPA
 * eigenvectors)
 * @param[in] descZ descriptor for Z
 * @param[in] jbegin first owned state, 1-based (checked against the
 * ownership derived internally from the 1D grid/NC; must match by
 * construction)
 * @param[in] ncols number of states owned by this rank (checked against
 * the internally derived ownership)
 * @param[in] panel_ld leading dimension of panel (= xNsize)
 * @param[in, out] panel this rank's 1D state-column panel buffer; on
 * return, panel[:, 0..ncols-1] (column-major, ld = panel_ld) holds the
 * full eigenvectors for this rank's owned 1-based states
 * [jbegin, jbegin+ncols-1]
 * @return 0 on success, -1 on failure (same value on all ranks)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int RedistBlockCyclicToStatePanel(long int xNsize,
                                  double complex *Z, int *descZ,
                                  long int jbegin, long int ncols,
                                  long int panel_ld, double complex *panel) {
  int i_negone = -1, i_zero_i = 0, info;
  int ictxt_1d, nprow_1, npcol_1, myrow_1, mycol_1;
  int desc1d[9];
  int lld;
  long int NC, mb1, nb1;
  const long int i_one = 1;
  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  NC = (xNsize + size - 1) / size;

  blacs_get_(&i_negone, &i_zero_i, &ictxt_1d);
  nprow_1 = 1; npcol_1 = size;
  blacs_gridinit_(&ictxt_1d, "R", &nprow_1, &npcol_1);
  blacs_gridinfo_(&ictxt_1d, &nprow_1, &npcol_1, &myrow_1, &mycol_1);

  /* Same ownership-consistency guard as RedistPanelToBlockCyclic: the
     1x P 'R' grid must map mycol to this rank 1:1, and the caller-
     supplied (jbegin, ncols) must agree with the ownership derived
     internally from mycol_1/NC/xNsize. Synchronize the verdict so no
     rank enters the collective pzgemr2d_ alone. */
  {
    int ok = (mycol_1 == myrank) ? 0 : -1, gok;
    long int nc_expect = (mycol_1 < (int)((xNsize + NC - 1) / NC))
                           ? (((long int)mycol_1 + 1) * NC <= xNsize
                                ? NC : xNsize - (long int)mycol_1 * NC)
                           : 0;
    long int jb_expect = (long int)mycol_1 * NC + 1;
    if (ok != 0) {
      fprintf(stdout,
              "  Error: BLACS 1D grid column (%d) does not match MPI rank (%d):\n"
              "         panel ownership is inconsistent; aborting redistribution.\n",
              mycol_1, myrank);
    }
    if (ncols != nc_expect ||
        (ncols > 0 && jbegin != jb_expect)) {
      ok = -1;
      fprintf(stdout,
              "  Error: caller-supplied panel ownership (jbegin=%ld, ncols=%ld) does\n"
              "         not match the internally derived ownership (jbegin=%ld, ncols=%ld)\n"
              "         for rank %d: panel ownership is inconsistent; aborting redistribution.\n",
              jbegin, ncols, jb_expect, nc_expect, myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      blacs_gridexit_(&ictxt_1d);
      return -1;
    }
  }

  lld = (panel_ld > 0) ? (int)panel_ld : 1;
  mb1 = xNsize;
  nb1 = NC;
  descinit_(desc1d, &xNsize, &xNsize, &mb1, &nb1, &i_zero_i, &i_zero_i,
            &ictxt_1d, &lld, &info);
  /* Synchronize the descinit_ verdict so no rank enters the collective
     pzgemr2d_ with an invalid descriptor while others proceed. */
  {
    int ok = (info == 0) ? 0 : -1, gok;
    if (ok != 0) {
      fprintf(stdout,
              "  Error: descinit_ failed (info=%d) for the state-panel descriptor\n"
              "         on rank %d; aborting redistribution.\n", info, myrank);
    }
    MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gok != 0) {
      blacs_gridexit_(&ictxt_1d);
      return -1;
    }
  }

  pzgemr2d_(&xNsize, &xNsize,
           Z, (long int *)&i_one, (long int *)&i_one, descZ,
           panel, (long int *)&i_one, (long int *)&i_one, desc1d,
           &descZ[1]);

  blacs_gridexit_(&ictxt_1d);
  return 0;
}

#endif

/* HPhi  -  Quantum Lattice Model Simulator */
#include "green_output.h"
#include "DefCommon.h"
#include "FileIO.h"
#include "global.h"
#include "wrapperMPI.h"
#include <string.h>
#include <stdlib.h>
#ifdef MPI
#include <mpi.h>
#endif

static int GreenOutputAggregateFamily(const struct BindStruct *X)
{
  if (X->Def.iOutputGreenFormat != OUTPUTGREENFORMAT_AGGREGATE) return 0;
  switch (X->Def.iCalcType) {
  case TPQCalc:
  case cTPQ:
    return 1;
  case TimeEvolution:
    return 2;
  case FullDiag:
  case CG:
    return 3;
  default:
    return 0;
  }
}

int GreenOutputUsesAggregate(const struct BindStruct *X)
{
  return GreenOutputAggregateFamily(X) != 0;
}

int GreenOutputKindUsesAggregate(const struct BindStruct *X, GreenOutputKind kind)
{
  if (!GreenOutputUsesAggregate(X)) return FALSE;
  return TRUE;
}

const char *GreenOutputOpenMode(const struct BindStruct *X)
{
  return GreenOutputUsesAggregate(X) ? "a" : "w";
}

static const char *GreenOutputFormatForFamily(int family, GreenOutputKind kind)
{
  if (family == 1) {
    switch (kind) {
    case GreenOutputOneBody: return cFileName1BGreen_TPQ_Aggregate;
    case GreenOutputTwoBody: return cFileName2BGreen_TPQ_Aggregate;
    case GreenOutputThreeBody: return cFileName3BGreen_TPQ_Aggregate;
    case GreenOutputFourBody: return cFileName4BGreen_TPQ_Aggregate;
    case GreenOutputSixBody: return cFileName6BGreen_TPQ_Aggregate;
    case GreenOutputNBody: return cFileNameNBodyG_TPQ_Aggregate;
    case GreenOutputAnomalous: return cFileNameAnomalousG_TPQ_Aggregate;
    }
  }
  if (family == 2) {
    switch (kind) {
    case GreenOutputOneBody: return cFileName1BGreen_TE_Aggregate;
    case GreenOutputTwoBody: return cFileName2BGreen_TE_Aggregate;
    case GreenOutputThreeBody: return cFileName3BGreen_TE_Aggregate;
    case GreenOutputFourBody: return cFileName4BGreen_TE_Aggregate;
    case GreenOutputSixBody: return cFileName6BGreen_TE_Aggregate;
    case GreenOutputNBody: return cFileNameNBodyG_TE_Aggregate;
    case GreenOutputAnomalous: return cFileNameAnomalousG_TE_Aggregate;
    }
  }
  if (family == 3) {
    switch (kind) {
    case GreenOutputOneBody: return cFileName1BGreen_Eigen_Aggregate;
    case GreenOutputTwoBody: return cFileName2BGreen_Eigen_Aggregate;
    case GreenOutputThreeBody: return cFileName3BGreen_Eigen_Aggregate;
    case GreenOutputFourBody: return cFileName4BGreen_Eigen_Aggregate;
    case GreenOutputSixBody: return cFileName6BGreen_Eigen_Aggregate;
    case GreenOutputNBody: return cFileNameNBodyG_Eigen_Aggregate;
    case GreenOutputAnomalous: return cFileNameAnomalousG_Eigen_Aggregate;
    }
  }
  return NULL;
}

int GreenOutputFileName(const struct BindStruct *X, GreenOutputKind kind, char *sdt)
{
  int family;
  const char *fmt = NULL;
  if (!GreenOutputKindUsesAggregate(X, kind)) return -1;
  family = GreenOutputAggregateFamily(X);
  fmt = GreenOutputFormatForFamily(family, kind);
  if (fmt == NULL) return -1;
  sprintf(sdt, fmt, X->Def.CDataFileHead);
  return 0;
}

int GreenOutputWriteIndexPrefix(FILE *fp, const struct BindStruct *X)
{
  switch (GreenOutputAggregateFamily(X)) {
  case 1:
    fprintf(fp, " %4d %4d ", X->Def.irand, X->Def.istep);
    return 0;
  case 2:
    fprintf(fp, " %4d ", X->Def.istep);
    return 0;
  case 3:
    fprintf(fp, " %4d ", X->Phys.eigen_num);
    return 0;
  default:
    return 0;
  }
}

static int GreenOutputTruncateFile(struct BindStruct *X, GreenOutputKind kind)
{
  FILE *fp = NULL;
  char sdt[D_FileNameMax];
  if (!GreenOutputKindUsesAggregate(X, kind)) return 0;
  if (GreenOutputFileName(X, kind, sdt) != 0) return -1;
  if (childfopenMPI(sdt, "w", &fp) != 0) return -1;
  fclose(fp);
  return 0;
}

int GreenOutputInitializeAggregateFiles(struct BindStruct *X)
{
  if (!GreenOutputUsesAggregate(X)) return 0;
  if (X->Def.NCisAjt > 0 &&
      GreenOutputTruncateFile(X, GreenOutputOneBody) != 0) return -1;
  if (X->Def.NCisAjtCkuAlvDC > 0 &&
      GreenOutputTruncateFile(X, GreenOutputTwoBody) != 0) return -1;
  if (X->Def.NTBody > 0 &&
      GreenOutputTruncateFile(X, GreenOutputThreeBody) != 0) return -1;
  if (X->Def.NFBody > 0 &&
      GreenOutputTruncateFile(X, GreenOutputFourBody) != 0) return -1;
  if (X->Def.NSBody > 0 &&
      GreenOutputTruncateFile(X, GreenOutputSixBody) != 0) return -1;
  if (X->Def.NNBodyG > 0 &&
      GreenOutputTruncateFile(X, GreenOutputNBody) != 0) return -1;
  if (X->Def.NAnomalousG > 0 &&
      GreenOutputTruncateFile(X, GreenOutputAnomalous) != 0) return -1;
  return 0;
}

int GreenOutputUsesTPQDataAggregate(const struct BindStruct *X)
{
  if (X->Def.iOutputGreenFormat != OUTPUTGREENFORMAT_AGGREGATE) return 0;
  return (X->Def.iCalcType == TPQCalc || X->Def.iCalcType == cTPQ);
}

static const char *GreenOutputTPQDataBaseName(GreenOutputTPQDataKind kind)
{
  switch (kind) {
  case GreenOutputTPQDataSS: return cFileNameSS_TPQ_Aggregate;
  case GreenOutputTPQDataNorm: return cFileNameNorm_TPQ_Aggregate;
  case GreenOutputTPQDataFlct: return cFileNameFlct_TPQ_Aggregate;
  }
  return NULL;
}

int GreenOutputTPQDataFileName(const struct BindStruct *X, GreenOutputTPQDataKind kind, char *sdt)
{
  const char *base = NULL;
  if (!GreenOutputUsesTPQDataAggregate(X)) return -1;
  base = GreenOutputTPQDataBaseName(kind);
  if (base == NULL) return -1;
  if (X->Def.iOutputDataHead == 1) {
    sprintf(sdt, "%s_%s", X->Def.CDataFileHead, base);
  } else {
    sprintf(sdt, "%s", base);
  }
  return 0;
}

static const char *GreenOutputTPQDataHeader(GreenOutputTPQDataKind kind)
{
  switch (kind) {
  case GreenOutputTPQDataSS:
    return " # set, step, inv_tmp, energy, phys_var, phys_doublon, phys_num\n";
  case GreenOutputTPQDataNorm:
    return " # set, step, inv_temp, global_norm, global_1st_norm\n";
  case GreenOutputTPQDataFlct:
    return " # set, step, inv_temp, N, N^2, D, D^2, Sz, Sz^2\n";
  }
  return NULL;
}

static int GreenOutputInitializeTPQDataFile(struct BindStruct *X, GreenOutputTPQDataKind kind)
{
  FILE *fp = NULL;
  char sdt[D_FileNameMax];
  const char *header = NULL;
  if (!GreenOutputUsesTPQDataAggregate(X)) return 0;
  if (GreenOutputTPQDataFileName(X, kind, sdt) != 0) return -1;
  header = GreenOutputTPQDataHeader(kind);
  if (header == NULL) return -1;
  if (childfopenMPI(sdt, "w", &fp) != 0) return -1;
  fprintf(fp, "%s", header);
  fclose(fp);
  return 0;
}

int GreenOutputInitializeTPQDataAggregateFiles(struct BindStruct *X)
{
  if (!GreenOutputUsesTPQDataAggregate(X)) return 0;
  if (GreenOutputInitializeTPQDataFile(X, GreenOutputTPQDataSS) != 0) return -1;
  if (GreenOutputInitializeTPQDataFile(X, GreenOutputTPQDataNorm) != 0) return -1;
  if (GreenOutputInitializeTPQDataFile(X, GreenOutputTPQDataFlct) != 0) return -1;
  return 0;
}

void GreenOutputWriteTPQSSRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp)
{
  if (GreenOutputUsesTPQDataAggregate(X)) {
    fprintf(fp, " %4d %4d %.16lf  %.16lf %.16lf %.16lf %.16lf\n",
            X->Def.irand, step, inv_temp, X->Phys.energy, X->Phys.var,
            X->Phys.doublon, X->Phys.num);
  } else {
    fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n",
            inv_temp, X->Phys.energy, X->Phys.var, X->Phys.doublon,
            X->Phys.num, step);
  }
}

void GreenOutputWriteTPQNormRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp,
                                double norm, double first_norm)
{
  if (GreenOutputUsesTPQDataAggregate(X)) {
    fprintf(fp, " %4d %4d %.16lf %.16lf %.16lf\n",
            X->Def.irand, step, inv_temp, norm, first_norm);
  } else {
    fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, norm, first_norm, step);
  }
}

void GreenOutputWriteTPQFlctRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp)
{
  if (GreenOutputUsesTPQDataAggregate(X)) {
    fprintf(fp, " %4d %4d %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n",
            X->Def.irand, step, inv_temp, X->Phys.num, X->Phys.num2,
            X->Phys.doublon, X->Phys.doublon2, X->Phys.Sz, X->Phys.Sz2);
  } else {
    fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n",
            inv_temp, X->Phys.num, X->Phys.num2, X->Phys.doublon,
            X->Phys.doublon2, X->Phys.Sz, X->Phys.Sz2, step);
  }
}

/* ------------------------------------------------------------------------
 * Partial-output / manifest-based aggregate merge (phase 3a, Mode 1).
 * See docs/superpowers/specs/2026-07-11-elpa-fulldiag-phase3-design.md §3
 * ("集約 Green ファイルのパーシャル出力") for the full design rationale.
 * ------------------------------------------------------------------------ */

/* Number of GreenOutputKind values tracked by the manifest. */
#define GREEN_OUTPUT_NKIND (GreenOutputAnomalous + 1)

/*
 * Manifest record for one (rank, kind) pair. Field order and widths are a
 * FIXED wire shape: int attempted, opened, open_error, closed_ok; long int
 * bytes; char part_path[256]; char final_path[256] -- because
 * GreenOutputMergePartials() transfers arrays of this struct between ranks
 * with MPI_Gather(..., sizeof(record), MPI_BYTE, ...). That is an in-memory,
 * same-binary, homogeneous-ABI transfer (every HPhi MPI rank in a run
 * executes the identical binary), NOT a portable/persisted wire format --
 * struct padding, endianness, and long-int width are whatever this
 * compilation produced, and this struct must never be written to disk or
 * exchanged across HPhi binaries/architectures.
 */
typedef struct {
  int attempted;
  int opened;
  int open_error;
  int closed_ok;
  long int bytes;
  char part_path[256];
  char final_path[256];
} GreenOutputManifestRecord;

static GreenOutputManifestRecord g_greenOutputManifest[GREEN_OUTPUT_NKIND];
/* Sticky "a close ever failed for this kind" flag, kept out of the gathered
   record itself (only closed_ok, derived from this, is part of the wire
   shape); reset together with the manifest by GreenOutputSetPartialSuffix(). */
static int g_greenOutputCloseFailed[GREEN_OUTPUT_NKIND];
static int g_greenOutputPartialActive = 0;
static int g_greenOutputPartialRank = 0;

void GreenOutputSetPartialSuffix(int rank)
{
  memset(g_greenOutputManifest, 0, sizeof(g_greenOutputManifest));
  memset(g_greenOutputCloseFailed, 0, sizeof(g_greenOutputCloseFailed));
  g_greenOutputPartialRank = rank;
  g_greenOutputPartialActive = 1;
}

void GreenOutputClearPartialSuffix(void)
{
  /* Deliberately does NOT touch g_greenOutputManifest: GreenOutputMergePartials()
     must still be able to read this session's records after the session is
     closed. Only the next GreenOutputSetPartialSuffix() zeroes it. */
  g_greenOutputPartialActive = 0;
}

/* Join a childfopenMPI()-relative path with the output-folder prefix, the
   same way childfopenMPI()/FileIO.c does internally, so remove() targets
   the exact on-disk path childfopenMPI() would open. */
static void GreenOutputJoinOutputPath(const char *rel, char *out, size_t outsz)
{
  out[0] = '\0';
  strncat(out, cParentOutputFolder, outsz - 1);
  strncat(out, rel, outsz - 1 - strlen(out));
}

int GreenOutputOpenAggregate(struct BindStruct *X, GreenOutputKind kind, FILE **fp)
{
  char final_sdt[D_FileNameMax];
  GreenOutputManifestRecord *rec;
  int n;

  if (fp == NULL) return -1;
  if (!GreenOutputKindUsesAggregate(X, kind)) return -1;
  if (GreenOutputFileName(X, kind, final_sdt) != 0) return -1;

  if (!g_greenOutputPartialActive) {
    return childfopenMPI(final_sdt, GreenOutputOpenMode(X), fp);
  }
  if (kind < 0 || kind > GreenOutputAnomalous) return -1;

  rec = &g_greenOutputManifest[kind];
  rec->attempted = 1;
  if (rec->open_error) return -1; /* sticky: this kind already failed this session */

  n = snprintf(rec->final_path, sizeof(rec->final_path), "%s", final_sdt);
  if (n < 0 || (size_t)n >= sizeof(rec->final_path)) { rec->open_error = 1; return -1; }
  n = snprintf(rec->part_path, sizeof(rec->part_path), "%s.part%d", final_sdt, g_greenOutputPartialRank);
  if (n < 0 || (size_t)n >= sizeof(rec->part_path)) { rec->open_error = 1; return -1; }

  if (!rec->opened) {
    /* First open of this kind in the session: unlink any stale part file
       (from a previous, possibly-failed run) before creating a fresh one. */
    char joined[sizeof(rec->part_path) + 64];
    GreenOutputJoinOutputPath(rec->part_path, joined, sizeof(joined));
    remove(joined); /* best-effort; ENOENT etc. are not errors here */
    if (childfopenMPI(rec->part_path, "w", fp) != 0) { rec->open_error = 1; return -1; }
    rec->opened = 1;
  } else {
    /* Later opens of the same kind in this session append. */
    if (childfopenMPI(rec->part_path, "a", fp) != 0) { rec->open_error = 1; return -1; }
  }
  return 0;
}

int GreenOutputCloseAggregate(GreenOutputKind kind, FILE *fp)
{
  if (fp == NULL) return -1;
  if (!g_greenOutputPartialActive) {
    return (fclose(fp) == 0) ? 0 : -1;
  }
  if (kind < 0 || kind > GreenOutputAnomalous) {
    fclose(fp);
    return -1;
  }
  {
    GreenOutputManifestRecord *rec = &g_greenOutputManifest[kind];
    long int pos = ftell(fp);
    int close_rc = fclose(fp);
    int ok = (pos >= 0 && close_rc == 0);
    if (ok) rec->bytes = pos;
    if (!ok) g_greenOutputCloseFailed[kind] = 1;
    rec->closed_ok = g_greenOutputCloseFailed[kind] ? 0 : 1;
    return ok ? 0 : -1;
  }
}

#ifdef MPI
int GreenOutputMergePartials(struct BindStruct *X)
{
  int nprocs_l = 1, myrank_l = 0;
  int rc = 0;
  int r, k;
  GreenOutputManifestRecord *all = NULL;

  (void)X; /* not currently needed: every path/kind is already resolved in the manifest */

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_l);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_l);

  if (myrank_l == 0) {
    all = (GreenOutputManifestRecord *)malloc(
        (size_t)nprocs_l * GREEN_OUTPUT_NKIND * sizeof(GreenOutputManifestRecord));
    if (all == NULL) {
      fprintf(stdoutMPI, "Error: GreenOutputMergePartials: manifest gather buffer allocation failed.\n");
      exitMPI(1); /* MPI_Abort()s the whole communicator; safe from a single rank */
    }
  }

  /* See the GreenOutputManifestRecord comment above: this is a same-binary,
     in-memory MPI_BYTE transfer, not a portable wire format. */
  MPI_Gather(g_greenOutputManifest,
             GREEN_OUTPUT_NKIND * (int)sizeof(GreenOutputManifestRecord), MPI_BYTE,
             all,
             GREEN_OUTPUT_NKIND * (int)sizeof(GreenOutputManifestRecord), MPI_BYTE,
             0, MPI_COMM_WORLD);

  if (myrank_l == 0) {
    /* Pass 1: an attempted record is trusted only if its whole lifecycle
       succeeded: the open did not fail, the file was actually opened, and
       the last close (fflush included) succeeded. A part that was opened
       but never closed, or whose fclose failed, may be missing buffered
       data even though the file itself is still readable. Never infer
       success from file sizes alone -- the manifest decides (pass 2 only
       cross-checks reality against the manifest's own claims). */
    for (r = 0; r < nprocs_l && rc == 0; r++) {
      for (k = 0; k < GREEN_OUTPUT_NKIND; k++) {
        GreenOutputManifestRecord *rr = &all[r * GREEN_OUTPUT_NKIND + k];
        if (rr->attempted &&
            (rr->open_error || !rr->opened || !rr->closed_ok)) { rc = -1; break; }
      }
    }
    /* Pass 2: verify every part file the manifest claims succeeded can
       still actually be opened for reading, before publishing anything.
       Catches a part file that vanished/was corrupted after a successful
       close, rather than discovering it mid-concatenation (which could
       otherwise leave a partially-published final file). */
    if (rc == 0) {
      for (r = 0; r < nprocs_l && rc == 0; r++) {
        for (k = 0; k < GREEN_OUTPUT_NKIND; k++) {
          GreenOutputManifestRecord *rr = &all[r * GREEN_OUTPUT_NKIND + k];
          FILE *probe = NULL;
          if (!rr->attempted) continue; /* legitimately empty / zero-owner rank */
          if (childfopenMPI(rr->part_path, "rb", &probe) != 0) { rc = -1; break; }
          /* Cross-check the on-disk length against the manifest's own
             recorded byte count (ftell at the last successful close): a
             mismatch means the part changed after the writer closed it. */
          if (fseek(probe, 0L, SEEK_END) != 0 || ftell(probe) != rr->bytes) {
            fclose(probe); rc = -1; break;
          }
          fclose(probe);
        }
      }
    }
  }

  if (myrank_l == 0 && rc == 0) {
    /* Two-phase publish (write-to-temp, then rename):
       Phase A concatenates each attempted kind's parts into a private
       <final>.tmp_merge file, with every fread/fwrite/ferror/fclose checked.
       Phase B, entered only if EVERY kind's temp was written and closed
       cleanly, rename()s each temp onto its final name. This way a write
       failure (disk full, quota, I/O error) during concatenation can never
       leave a truncated file under the final name, and no final file is
       published unless all kinds succeeded. Renaming per kind is the
       Mode-1 replacement for GreenOutputInitializeAggregateFiles()
       (which Mode 1 must not call directly). */
    char tmp_joined[GREEN_OUTPUT_NKIND][sizeof(((GreenOutputManifestRecord *)0)->final_path) + 80];
    char final_joined[GREEN_OUTPUT_NKIND][sizeof(((GreenOutputManifestRecord *)0)->final_path) + 64];
    int kind_active[GREEN_OUTPUT_NKIND];

    for (k = 0; k < GREEN_OUTPUT_NKIND; k++) kind_active[k] = 0;

    /* Phase A: concatenate into temp files. */
    for (k = 0; k < GREEN_OUTPUT_NKIND && rc == 0; k++) {
      int any_attempted = 0;
      int n;
      FILE *fout = NULL;

      for (r = 0; r < nprocs_l; r++) {
        GreenOutputManifestRecord *rr = &all[r * GREEN_OUTPUT_NKIND + k];
        if (rr->attempted) {
          any_attempted = 1;
          GreenOutputJoinOutputPath(rr->final_path, final_joined[k], sizeof(final_joined[k]));
          break;
        }
      }
      if (!any_attempted) continue;

      n = snprintf(tmp_joined[k], sizeof(tmp_joined[k]), "%s.tmp_merge", final_joined[k]);
      if (n < 0 || (size_t)n >= sizeof(tmp_joined[k])) { rc = -1; break; }
      /* Only mark active once tmp_joined[k] holds a valid path -- the
         failure-cleanup path remove()s every active kind's temp path. */
      kind_active[k] = 1;

      fout = fopen(tmp_joined[k], "wb");
      if (fout == NULL) { rc = -1; break; }

      for (r = 0; r < nprocs_l && rc == 0; r++) {
        GreenOutputManifestRecord *rr = &all[r * GREEN_OUTPUT_NKIND + k];
        FILE *fin = NULL;
        char part_joined[sizeof(rr->part_path) + 64];
        char buf[8192];
        size_t got;

        if (!rr->attempted) continue; /* zero-owner rank: contributes nothing */
        GreenOutputJoinOutputPath(rr->part_path, part_joined, sizeof(part_joined));
        fin = fopen(part_joined, "rb");
        if (fin == NULL) { rc = -1; break; }
        while ((got = fread(buf, 1, sizeof(buf), fin)) > 0) {
          if (fwrite(buf, 1, got, fout) != got) { rc = -1; break; }
        }
        if (rc == 0 && ferror(fin)) rc = -1; /* short read due to I/O error */
        if (fclose(fin) != 0) rc = -1;
      }
      if (fclose(fout) != 0) rc = -1; /* flush failure = truncated temp */
    }

    /* Phase B: publish by rename, only if every kind concatenated cleanly. */
    if (rc == 0) {
      for (k = 0; k < GREEN_OUTPUT_NKIND && rc == 0; k++) {
        if (!kind_active[k]) continue;
        if (rename(tmp_joined[k], final_joined[k]) != 0) rc = -1;
      }
    }

    if (rc == 0) {
      /* Only ever delete part files once every kind published successfully. */
      for (r = 0; r < nprocs_l; r++) {
        for (k = 0; k < GREEN_OUTPUT_NKIND; k++) {
          GreenOutputManifestRecord *rr = &all[r * GREEN_OUTPUT_NKIND + k];
          if (rr->attempted) {
            char joined[sizeof(rr->part_path) + 64];
            GreenOutputJoinOutputPath(rr->part_path, joined, sizeof(joined));
            remove(joined);
          }
        }
      }
    } else {
      /* Failure: remove whatever temps exist (best effort), keep all part
         files on disk for diagnosis, publish nothing further. (If a rename
         in Phase B failed partway, kinds already renamed stay published --
         rename is the smallest possible window -- but rc=-1 is still
         returned everywhere and no part file is deleted.) */
      for (k = 0; k < GREEN_OUTPUT_NKIND; k++) {
        if (kind_active[k]) remove(tmp_joined[k]);
      }
      fprintf(stdoutMPI,
              "Error: GreenOutputMergePartials: merge failed; partial (.part*) files are kept for diagnosis.\n");
    }
  }

  /* Single rendezvous, LAST: every rank reaches this unconditionally (no
     early return exists between the Gather above and here on any rank), so
     the verdict broadcast from rank 0 covers manifest-validation failures
     AND publish-phase (Phase A/B) failures alike -- on any failure, every
     rank returns nonzero. */
  MPI_Bcast(&rc, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (myrank_l == 0 && all != NULL) free(all);
  return rc;
}
#else
int GreenOutputMergePartials(struct BindStruct *X)
{
  (void)X;
  return 0; /* no MPI: no cross-rank partial files exist to merge; no-op */
}
#endif

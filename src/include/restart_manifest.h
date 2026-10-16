#ifndef HPHI_RESTART_MANIFEST_H
#define HPHI_RESTART_MANIFEST_H

#include <stdint.h>
#include <stdio.h>

struct BindStruct;

struct RestartManifestRow {
  uint64_t set, ranks, gen, step, ncomp;
};

struct RestartManifest {
  int method;
  int basis;
  size_t nrows, capacity;
  struct RestartManifestRow *rows;
};

int RestartManifestParse(FILE *fp, struct RestartManifest *m);
int RestartManifestFormat(const struct RestartManifest *m, FILE *fp);
const struct RestartManifestRow *RestartManifestFind(
    const struct RestartManifest *m, uint64_t set);
int RestartManifestUpsert(
    struct RestartManifest *m, const struct RestartManifestRow *row);
void RestartManifestFree(struct RestartManifest *m);

/* Pure collision-free generation selection; fails on mixed published parity. */
int RestartChooseGeneration(
    const struct RestartManifest *old, const uint64_t *sets, int nset,
    uint64_t *gen);

int RestartManifestLoad(
    const struct BindStruct *X, struct RestartManifest *m, int *exists);
int RestartManifestPublish(
    const struct BindStruct *X, const struct RestartManifest *m);

#endif

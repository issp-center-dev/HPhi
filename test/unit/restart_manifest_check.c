#include <stdio.h>
#include <string.h>

#include "restart_manifest.h"
#include "DefCommon.h"
#include "restart_check.h"

static FILE *text_file(const char *text)
{
  FILE *fp = tmpfile();
  if (fp == NULL) return NULL;
  if (fputs(text, fp) == EOF) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  return fp;
}

static int check_round_trip(void)
{
  struct RestartManifest m, again;
  FILE *fp = text_file(
      "hphi_restart_manifest 1\n"
      "method TPQ\n"
      "basis raw\n"
      "set 0 ranks 4 gen 1 step 20 ncomp 1\n"
      "set 1 ranks 2 gen 0 step 13 ncomp 1\n");
  CHECK(fp != NULL);
  int status = RestartManifestParse(fp, &m);
  CHECK(fclose(fp) == 0);
  CHECK(status == 0);
  CHECK(m.method == TPQCalc && m.basis == 0 && m.nrows == 2);
  CHECK(RestartManifestFind(&m, 1) != NULL);
  CHECK(RestartManifestFind(&m, 1)->ranks == 2);
  CHECK(RestartManifestFind(&m, 1)->step == 13);
  CHECK(RestartManifestFind(&m, 7) == NULL);

  struct RestartManifestRow row = {1, 2, 1, 20, 1};
  CHECK(RestartManifestUpsert(&m, &row) == 0);
  CHECK(m.nrows == 2 && RestartManifestFind(&m, 1)->gen == 1);
  row.set = 2;
  CHECK(RestartManifestUpsert(&m, &row) == 0);
  CHECK(m.nrows == 3);

  fp = tmpfile();
  CHECK(fp != NULL);
  CHECK(RestartManifestFormat(&m, fp) == 0);
  rewind(fp);
  status = RestartManifestParse(fp, &again);
  CHECK(fclose(fp) == 0);
  CHECK(status == 0);
  CHECK(again.nrows == 3 && RestartManifestFind(&again, 2)->ranks == 2);
  RestartManifestFree(&m);
  RestartManifestFree(&again);
  return 0;
}

static int check_names(void)
{
  const int methods[] = {CG, TPQCalc, cTPQ, TimeEvolution};
  const char *names[] = {"CG", "TPQ", "cTPQ", "TimeEvolution"};
  for (int i = 0; i < 4; ++i) {
    char text[128];
    snprintf(text, sizeof(text),
             "hphi_restart_manifest 1\nmethod %s\nbasis sector\n", names[i]);
    FILE *fp = text_file(text);
    struct RestartManifest m;
    CHECK(fp != NULL);
    CHECK(RestartManifestParse(fp, &m) == 0);
    CHECK(fclose(fp) == 0);
    CHECK(m.method == methods[i] && m.basis == 1);
    RestartManifestFree(&m);
  }
  return 0;
}

static int check_bad_inputs(void)
{
  const char *bad[] = {
    "hphi_restart_manifest 2\nmethod CG\nbasis raw\n",
    "hphi_restart_manifest 1\nmethod CG\nbasis raw\nbogus 1\n",
    "hphi_restart_manifest 1\nmethod CG\nbasis raw\nset 0 ranks 1 gen 2 step 1 ncomp 2\n",
    "hphi_restart_manifest 1\nmethod CG\nbasis raw\nset 0 ranks 1 gen 0 step 1 ncomp 2\nset 0 ranks 1 gen 0 step 1 ncomp 2\n",
    "hphi_restart_manifest 1\nmethod Lanczos\nbasis raw\n",
    "hphi_restart_manifest 1\nmethod CG\nbasis other\n",
    "hphi_restart_manifest 1\nmethod CG\nbasis raw\nset 0 ranks 0 gen 0 step 1 ncomp 1\n",
    "hphi_restart_manifest 1\nmethod CG\nbasis raw\nset 0 ranks 1 gen 0 step 1 ncomp 0\n"
  };
  for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
    FILE *fp = text_file(bad[i]);
    struct RestartManifest m;
    CHECK(fp != NULL);
    int status = RestartManifestParse(fp, &m);
    CHECK(fclose(fp) == 0);
    CHECK(status != 0);
    RestartManifestFree(&m);
  }
  return 0;
}

int main(void)
{
  CHECK(check_round_trip() == 0);
  CHECK(check_names() == 0);
  CHECK(check_bad_inputs() == 0);
  /* Add to restart_manifest_check.c main in Task 9 (V1/V7); free is local to this block. */
  {
    struct RestartManifest old = {0}; uint64_t sets[2] = {0, 1}, gen = 99;
    struct RestartManifestRow a = {0, 4, 0, 8, 1}, b = {1, 4, 1, 8, 1};
    CHECK(RestartChooseGeneration(&old, sets, 2, &gen) == 0 && gen == 0);
    CHECK(RestartManifestUpsert(&old, &a) == 0);
    CHECK(RestartChooseGeneration(&old, sets, 2, &gen) == 0 && gen == 1);
    CHECK(RestartManifestUpsert(&old, &b) == 0);
    CHECK(RestartChooseGeneration(&old, sets, 2, &gen) == -1); /* mixed parity */
    CHECK(RestartChooseGeneration(&old, sets, 1, &gen) == 0 && gen == 1);
    RestartManifestFree(&old);
  }

  puts("restart_manifest_check: ok");
  return 0;
}

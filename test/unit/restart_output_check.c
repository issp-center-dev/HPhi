#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "restart_output.h"
#include "restart_check.h"

static int write_file(const char *path, const char *text)
{
  FILE *fp = fopen(path, "w");
  if (!fp) return -1;
  fputs(text, fp);
  return fclose(fp);
}

static char *read_file(const char *path, char *buf, size_t n)
{
  FILE *fp = fopen(path, "r");
  size_t k = fp ? fread(buf, 1, n - 1, fp) : 0;
  buf[k] = 0;
  if (fp) fclose(fp);
  return buf;
}

int main(void)
{
  static char buf[1 << 16];
  char msg[256];
  int status;

  CHECK(write_file("cg.dat", "    Step   Residual-2-norm     Threshold      Energy\n"
                           "        1  1e-1 1e-5 -3.0\n"
                           "        2  1e-2 1e-5 -3.1\n"
                           "        3  1e-3 1e-5 -3.2\n"
                           "        4  1e-4") == 0);
  struct RestartSeriesRule cg_rule = {0, 3};
  struct RestartSeriesSpec cg = {"Step", -1, 0, 1, &cg_rule, 0};
  status = RestartRecoverTimeSeriesFile("cg.dat", &cg, msg, sizeof(msg));
  CHECK(status == 0);
  CHECK(strcmp(read_file("cg.dat", buf, sizeof(buf)),
             "    Step   Residual-2-norm     Threshold      Energy\n"
             "        1  1e-1 1e-5 -3.0\n"
             "        2  1e-2 1e-5 -3.1\n") == 0);

  /* Valid rows with 64 energy/Ritz values exceed 1023 bytes and survive. */
  {
    FILE *fp = fopen("long.dat", "w");
    CHECK(fp != NULL);
    fputs("    Step   Residual-2-norm     Threshold      Energy\n", fp);
    for (int row = 1; row <= 2; ++row) {
      fprintf(fp, "%9d %15.5e %15.5e      ", row, 1e-3, 1e-5);
      for (int e = 0; e < 64; ++e) fprintf(fp, " %15.5e", -1.0 - e);
      fputs("\n", fp);
    }
    CHECK(fclose(fp) == 0);
  }
  struct RestartSeriesRule long_rule = {0, 5};
  struct RestartSeriesSpec lg = {"Step", -1, 0, 1, &long_rule, 0};
  status = RestartRecoverTimeSeriesFile("long.dat", &lg, msg, sizeof(msg));
  CHECK(status == 0);
  read_file("long.dat", buf, sizeof(buf));
  {
    int lines = 0;
    for (char *p = buf; *p; ++p) lines += *p == '\n';
    CHECK(lines == 3);
  }

  /* A comparably long corrupt interior row must preserve the original file. */
  {
    FILE *fp = fopen("long-bad.dat", "w");
    CHECK(fp != NULL);
    fputs(" # h\n0 0 0", fp);
    for (int i = 0; i < 100; ++i) fputs(" 12345678901234567890", fp);
    fputs("\ngarbage", fp);
    for (int i = 0; i < 100; ++i) fputs(" 12345678901234567890", fp);
    fputs("\n0 2 0", fp);
    for (int i = 0; i < 100; ++i) fputs(" 12345678901234567890", fp);
    fputs("\n", fp);
    CHECK(fclose(fp) == 0);
  }
  read_file("long-bad.dat", buf, sizeof(buf));
  char original[1 << 16];
  memcpy(original, buf, sizeof(original));
  struct RestartSeriesRule bad_rule = {0, 2};
  struct RestartSeriesSpec bad = {"#", -1, 0, 1, &bad_rule, 0};
  status = RestartRecoverTimeSeriesFile("long-bad.dat", &bad, msg, sizeof(msg));
  CHECK(status != 0 && strstr(msg, "line 3") != NULL);
  CHECK(strcmp(read_file("long-bad.dat", buf, sizeof(buf)), original) == 0);

  /* A long EOF fragment is dropped only after EOF confirms no newline. */
  {
    FILE *fp = fopen("long-tail.dat", "w");
    CHECK(fp != NULL);
    fputs(" # h\n0 0 0", fp);
    for (int i = 0; i < 200; ++i) fputs(" 12345678901234567890", fp);
    CHECK(fclose(fp) == 0);
  }
  status = RestartRecoverTimeSeriesFile("long-tail.dat", &bad, msg, sizeof(msg));
  CHECK(status == 0);
  CHECK(strcmp(read_file("long-tail.dat", buf, sizeof(buf)), " # h\n") == 0);

  CHECK(write_file("bad.dat", " # h\n0.0 -1 0 0 0 0\n"
                              "garbage line\n0.2 -1 0 0 0 2\n") == 0);
  struct RestartSeriesRule te_rule = {0, 2};
  struct RestartSeriesSpec te = {"#", -1, -1, 1, &te_rule, 0};
  struct RestartSeriesRule agg_rules[] = {{0, 4}, {1, 3}};
  struct RestartSeriesSpec agg = {"#", 0, 1, 2, agg_rules, 1};
  status = RestartRecoverTimeSeriesFile("bad.dat", &te, msg, sizeof(msg));
  CHECK(status != 0 && strstr(msg, "line 3") != NULL);
  CHECK(strcmp(read_file("bad.dat", buf, sizeof(buf)),
             " # h\n0.0 -1 0 0 0 0\n"
             "garbage line\n0.2 -1 0 0 0 2\n") == 0);

  /* Signs are not unsigned decimal set/step values; each rejection preserves bytes. */
  CHECK(write_file("negative-step.dat", " # h\n0 -1\n") == 0);
  status = RestartRecoverTimeSeriesFile("negative-step.dat", &te, msg, sizeof(msg));
  CHECK(status != 0 && strstr(msg, "line 2") != NULL);
  CHECK(strcmp(read_file("negative-step.dat", buf, sizeof(buf)), " # h\n0 -1\n") == 0);

  CHECK(write_file("negative-set.dat", " # h\n-1 0\n") == 0);
  status = RestartRecoverTimeSeriesFile("negative-set.dat", &agg, msg, sizeof(msg));
  CHECK(status != 0 && strstr(msg, "line 2") != NULL);
  CHECK(strcmp(read_file("negative-set.dat", buf, sizeof(buf)), " # h\n-1 0\n") == 0);

  CHECK(write_file("blank-interior.dat", " # h\n0 0\n\n0 1\n") == 0);
  status = RestartRecoverTimeSeriesFile("blank-interior.dat", &te, msg, sizeof(msg));
  CHECK(status != 0 && strstr(msg, "line 3") != NULL);
  CHECK(strcmp(read_file("blank-interior.dat", buf, sizeof(buf)),
             " # h\n0 0\n\n0 1\n") == 0);

  CHECK(write_file("agg.dat", " # set, step, inv_tmp, energy\n"
                               " 0 0 0 1\n 1 0 0 1\n 0 3 1 1\n"
                               " 0 4 1 1\n 1 2 1 1\n 2 0 0 1\n") == 0);
  status = RestartRecoverTimeSeriesFile("agg.dat", &agg, msg, sizeof(msg));
  CHECK(status == 0);
  CHECK(strcmp(read_file("agg.dat", buf, sizeof(buf)),
             " # set, step, inv_tmp, energy\n 0 0 0 1\n 1 0 0 1\n"
             " 0 3 1 1\n 1 2 1 1\n") == 0);

  CHECK(write_file("ss.dat", " # time, energy, var, doublon, num, step\n"
                              "0.0 -1 0 0 0 0\n0.1 -1 0 0 0 1\n"
                              "0.2 -1 0 0 0 2\n") == 0);
  status = RestartRecoverTimeSeriesFile("ss.dat", &te, msg, sizeof(msg));
  CHECK(status == 0);
  CHECK(strcmp(read_file("ss.dat", buf, sizeof(buf)),
             " # time, energy, var, doublon, num, step\n"
             "0.0 -1 0 0 0 0\n0.1 -1 0 0 0 1\n") == 0);

  /* Failure to create the temporary file leaves the original untouched. */
  CHECK(write_file("temp-fail.dat", " # h\n0 0\n") == 0);
  CHECK(mkdir("temp-fail.dat.part", 0700) == 0);
  status = RestartRecoverTimeSeriesFile("temp-fail.dat", &te, msg, sizeof(msg));
  CHECK(status != 0);
  CHECK(strcmp(read_file("temp-fail.dat", buf, sizeof(buf)), " # h\n0 0\n") == 0);
  CHECK(rmdir("temp-fail.dat.part") == 0);

  /* Direct I/O decision points preserve the old bytes and remove .part. */
  const char *io_failures[] = {"open", "write", "close", "rename"};
  for (size_t i = 0; i < sizeof(io_failures)/sizeof(io_failures[0]); ++i) {
    CHECK(write_file("io-fail.dat", " # h\n0 0\n0 1\n") == 0);
    setenv("HPHI_TEST_RESTART_OUTPUT_FAIL", io_failures[i], 1);
    status = RestartRecoverTimeSeriesFile("io-fail.dat", &te, msg, sizeof(msg));
    unsetenv("HPHI_TEST_RESTART_OUTPUT_FAIL");
    CHECK(status != 0);
    CHECK(strcmp(read_file("io-fail.dat", buf, sizeof(buf)),
                 " # h\n0 0\n0 1\n") == 0);
    CHECK(access("io-fail.dat.part", F_OK) != 0);
  }

  remove("absent.dat");
  status = RestartRecoverTimeSeriesFile("absent.dat", &te, msg, sizeof(msg));
  CHECK(status == 0);
  puts("restart_output_check: ok");
  return 0;
}

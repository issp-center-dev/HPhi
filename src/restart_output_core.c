#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "restart_output.h"

#ifdef HPHI_TEST_RESTART_IO
static int output_test_fail(const char *stage)
{
  const char *wanted = getenv("HPHI_TEST_RESTART_OUTPUT_FAIL");
  return wanted != NULL && strcmp(wanted, stage) == 0;
}
#endif

/* Read one line of any length.  A line without a final newline is marked
 * incomplete so callers can discard an EOF fragment after seeing EOF. */
static int read_line(FILE *fp, char **buf, size_t *cap, size_t *len,
                     int *complete)
{
  *len = 0;
  *complete = 0;
  for (;;) {
    int ch = fgetc(fp);
    if (ch == EOF) {
      if (ferror(fp)) return -1;
      if (*len == 0) return 0;
      break;
    }
    if (*len + 2 > *cap) {
      size_t ncap = *cap ? 2 * *cap : 4096;
      char *nbuf = realloc(*buf, ncap);
      if (nbuf == NULL) return -1;
      *buf = nbuf;
      *cap = ncap;
    }
    (*buf)[(*len)++] = (char)ch;
    if (ch == '\n') {
      *complete = 1;
      break;
    }
  }
  if (*cap == 0) {
    *buf = malloc(1);
    if (*buf == NULL) return -1;
    *cap = 1;
  }
  (*buf)[*len] = '\0';
  return 1;
}

static int parse_u64(const char *token, size_t n, uint64_t *value)
{
  uint64_t parsed = 0;

  if (n == 0) return -1;
  for (size_t i = 0; i < n; ++i) {
    unsigned char ch = (unsigned char)token[i];
    uint64_t digit;
    if (ch < '0' || ch > '9') return -1;
    digit = (uint64_t)(ch - '0');
    if (parsed > (UINT64_MAX - digit) / 10) return -1;
    parsed = parsed * 10 + digit;
  }
  *value = parsed;
  return 0;
}

/* Returns 1 keep, 0 drop, -1 unparsable, 2 blank.  Tokens are scanned in
 * place so rows have no fixed maximum width. */
static int classify_data_row(const char *line,
                             const struct RestartSeriesSpec *spec)
{
  const char *tokens[2] = {NULL, NULL};
  size_t lengths[2] = {0, 0};
  const char *last = NULL;
  size_t last_len = 0;
  int ntokens = 0;

  for (const char *p = line; *p != '\0';) {
    while (*p != '\0' && isspace((unsigned char)*p)) ++p;
    if (*p == '\0') break;
    const char *start = p;
    while (*p != '\0' && !isspace((unsigned char)*p)) ++p;
    if (ntokens == spec->set_column) {
      tokens[0] = start;
      lengths[0] = (size_t)(p - start);
    }
    if (ntokens == spec->step_column) {
      tokens[1] = start;
      lengths[1] = (size_t)(p - start);
    }
    last = start;
    last_len = (size_t)(p - start);
    ++ntokens;
  }
  if (ntokens == 0) return 2;
  if (spec->step_column < 0) {
    tokens[1] = last;
    lengths[1] = last_len;
  }
  if (tokens[1] == NULL ||
      (spec->set_column >= 0 && tokens[0] == NULL)) return -1;

  uint64_t step, set = 0;
  if (parse_u64(tokens[1], lengths[1], &step) != 0) return -1;
  if (spec->set_column >= 0 &&
      parse_u64(tokens[0], lengths[0], &set) != 0) return -1;

  const struct RestartSeriesRule *rule = NULL;
  if (spec->set_column < 0) {
    if (spec->nrule > 0 && spec->rules != NULL) rule = &spec->rules[0];
  } else {
    for (int i = 0; i < spec->nrule; ++i) {
      if (spec->rules != NULL && spec->rules[i].set == set) {
        rule = &spec->rules[i];
        break;
      }
    }
  }
  if (rule == NULL) return spec->drop_unlisted_sets ? 0 : 1;
  return step < rule->keep_below_step ? 1 : 0;
}

static int is_header(const char *line, const struct RestartSeriesSpec *spec)
{
  while (*line == ' ' || *line == '\t') ++line;
  return spec->header_prefix != NULL &&
         strncmp(line, spec->header_prefix,
                 strlen(spec->header_prefix)) == 0;
}

static void set_message(char *message, size_t size, const char *format,
                        const char *path)
{
  if (message != NULL && size != 0) snprintf(message, size, format, path);
}

int RestartRecoverTimeSeriesFile(const char *path,
                                 const struct RestartSeriesSpec *spec,
                                 char *message, size_t size)
{
  FILE *in = NULL;
  FILE *out = NULL;
  char part[4096];
  char *buf = NULL;
  size_t cap = 0, len = 0;
  int complete = 0, lineno = 0, status = 0, got = 0, data_seen = 0;

  if (message != NULL && size != 0) message[0] = '\0';
  if (path == NULL || spec == NULL || spec->header_prefix == NULL ||
      spec->step_column < -1 || spec->set_column < -1 || spec->nrule < 0 ||
      (spec->nrule != 0 && spec->rules == NULL)) {
    set_message(message, size, "invalid output recovery arguments for %s",
                path != NULL ? path : "(null)");
    return -1;
  }

  errno = 0;
#ifdef HPHI_TEST_RESTART_IO
  if (output_test_fail("open")) {
    errno = EACCES;
    in = NULL;
  } else
#endif
  in = fopen(path, "rb");
  if (in == NULL) {
    if (errno == ENOENT) return 0;
    set_message(message, size, "cannot open %s", path);
    return -1;
  }
  if (snprintf(part, sizeof(part), "%s.part", path) >= (int)sizeof(part)) {
    fclose(in);
    set_message(message, size, "path too long for %s", path);
    return -1;
  }
  out = fopen(part, "wb");
  if (out == NULL) {
    fclose(in);
    set_message(message, size, "cannot create %s", part);
    return -1;
  }

  while ((got = read_line(in, &buf, &cap, &len, &complete)) == 1) {
    ++lineno;
    if (!complete) break;
    int keep;
    if (is_header(buf, spec)) {
      keep = 1;
    } else {
      int classification = classify_data_row(buf, spec);
      if (classification < 0 || (classification == 2 && data_seen)) {
        if (message != NULL && size != 0)
          snprintf(message, size, "unparsable data row at line %d of %s",
                   lineno, path);
        status = -1;
        break;
      }
      if (classification != 2) data_seen = 1;
      keep = classification != 0;
    }
    if (keep && (
#ifdef HPHI_TEST_RESTART_IO
        output_test_fail("write") ||
#endif
        fputs(buf, out) == EOF)) {
      set_message(message, size, "write failed for %s", part);
      status = -1;
      break;
    }
  }
  if (got < 0) {
    set_message(message, size, "read failed for %s", path);
    status = -1;
  }
  free(buf);
  if (fclose(in) != 0 && status == 0) {
    set_message(message, size, "close failed for %s", path);
    status = -1;
  }
  int close_failed = fclose(out) != 0;
#ifdef HPHI_TEST_RESTART_IO
  close_failed |= output_test_fail("close");
#endif
  if (close_failed && status == 0) {
    set_message(message, size, "close failed for %s", part);
    status = -1;
  }
  if (status != 0) {
    remove(part);
    return -1;
  }
#ifdef HPHI_TEST_RESTART_IO
  if (output_test_fail("rename") || rename(part, path) != 0) {
#else
  if (rename(part, path) != 0) {
#endif
    set_message(message, size, "cannot replace %s", path);
    remove(part);
    return -1;
  }
  return 0;
}

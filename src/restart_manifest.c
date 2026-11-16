#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "restart_manifest.h"
#include "DefCommon.h"

static const struct {
  const char *name;
  int type;
} methods[] = {
  {"CG", CG},
  {"TPQ", TPQCalc},
  {"cTPQ", cTPQ},
  {"TimeEvolution", TimeEvolution}
};

static const char *method_name(int type)
{
  for (size_t i = 0; i < sizeof(methods) / sizeof(methods[0]); ++i)
    if (methods[i].type == type) return methods[i].name;
  return NULL;
}

static int method_type(const char *name)
{
  for (size_t i = 0; i < sizeof(methods) / sizeof(methods[0]); ++i)
    if (strcmp(methods[i].name, name) == 0) return methods[i].type;
  return -1;
}

static int at_line_end(const char *line, int offset)
{
  const unsigned char *p = (const unsigned char *)line + offset;
  while (*p != '\0' && isspace(*p)) ++p;
  return *p == '\0';
}

static int parse_u64(const char *text, uint64_t *value)
{
  char *end;
  uintmax_t parsed;
  if (text == NULL || text[0] == '\0' || text[0] == '+' || text[0] == '-')
    return -1;
  errno = 0;
  parsed = strtoumax(text, &end, 10);
  if (errno == ERANGE || *end != '\0' || parsed > UINT64_MAX) return -1;
  *value = (uint64_t)parsed;
  return 0;
}

int RestartManifestParse(FILE *fp, struct RestartManifest *m)
{
  char line[256], key[32], value[64];
  int have_version = 0, have_method = 0, have_basis = 0;
  if (m == NULL) return -1;
  memset(m, 0, sizeof(*m));
  m->method = -1;
  m->basis = -1;
  if (fp == NULL) return -1;

  while (fgets(line, sizeof(line), fp) != NULL) {
    int offset = 0;
    if (strchr(line, '\n') == NULL && !feof(fp)) return -1;
    if (sscanf(line, "%31s", key) != 1) continue;
    if (strcmp(key, "hphi_restart_manifest") == 0) {
      uint64_t version;
      if (have_version || sscanf(line, "%*s %63s %n", value, &offset) != 1 ||
          !at_line_end(line, offset) || parse_u64(value, &version) != 0 ||
          version != 1)
        return -1;
      have_version = 1;
    } else if (strcmp(key, "method") == 0) {
      if (have_method || sscanf(line, "%*s %63s %n", value, &offset) != 1 ||
          !at_line_end(line, offset) || (m->method = method_type(value)) < 0)
        return -1;
      have_method = 1;
    } else if (strcmp(key, "basis") == 0) {
      if (have_basis || sscanf(line, "%*s %63s %n", value, &offset) != 1 ||
          !at_line_end(line, offset))
        return -1;
      if (strcmp(value, "raw") == 0)
        m->basis = 0;
      else if (strcmp(value, "sector") == 0)
        m->basis = 1;
      else
        return -1;
      have_basis = 1;
    } else if (strcmp(key, "set") == 0) {
      struct RestartManifestRow row;
      char fields[5][64];
      if (sscanf(line,
                 "set %63s ranks %63s gen %63s step %63s ncomp %63s %n",
                 fields[0], fields[1], fields[2], fields[3], fields[4],
                 &offset) != 5 ||
          !at_line_end(line, offset) ||
          parse_u64(fields[0], &row.set) != 0 ||
          parse_u64(fields[1], &row.ranks) != 0 ||
          parse_u64(fields[2], &row.gen) != 0 ||
          parse_u64(fields[3], &row.step) != 0 ||
          parse_u64(fields[4], &row.ncomp) != 0)
        return -1;
      if (row.ranks < 1 || row.gen > 1 || row.ncomp < 1 || row.ncomp > 2 ||
          RestartManifestFind(m, row.set) != NULL ||
          RestartManifestUpsert(m, &row) != 0)
        return -1;
    } else {
      return -1;
    }
  }
  if (ferror(fp)) return -1;
  return have_version && have_method && have_basis ? 0 : -1;
}

int RestartManifestFormat(const struct RestartManifest *m, FILE *fp)
{
  const char *name;
  if (m == NULL || fp == NULL) return -1;
  name = method_name(m->method);
  if (name == NULL || (m->basis != 0 && m->basis != 1) ||
      (m->nrows != 0 && m->rows == NULL))
    return -1;
  if (fprintf(fp, "hphi_restart_manifest 1\nmethod %s\nbasis %s\n", name,
              m->basis ? "sector" : "raw") < 0)
    return -1;
  for (size_t i = 0; i < m->nrows; ++i) {
    const struct RestartManifestRow *r = &m->rows[i];
    if (r->ranks < 1 || r->gen > 1 || r->ncomp < 1 || r->ncomp > 2)
      return -1;
    for (size_t j = 0; j < i; ++j)
      if (m->rows[j].set == r->set) return -1;
    if (fprintf(fp,
                "set %" PRIu64 " ranks %" PRIu64 " gen %" PRIu64
                " step %" PRIu64 " ncomp %" PRIu64 "\n",
                r->set, r->ranks, r->gen, r->step, r->ncomp) < 0)
      return -1;
  }
  return 0;
}

const struct RestartManifestRow *RestartManifestFind(
    const struct RestartManifest *m, uint64_t set)
{
  if (m == NULL) return NULL;
  for (size_t i = 0; i < m->nrows; ++i)
    if (m->rows[i].set == set) return &m->rows[i];
  return NULL;
}

int RestartManifestUpsert(
    struct RestartManifest *m, const struct RestartManifestRow *row)
{
  if (m == NULL || row == NULL) return -1;
  for (size_t i = 0; i < m->nrows; ++i) {
    if (m->rows[i].set == row->set) {
      m->rows[i] = *row;
      return 0;
    }
  }
  if (m->nrows == m->capacity) {
    size_t capacity;
    struct RestartManifestRow *rows;
    if (m->capacity > SIZE_MAX / 2) return -1;
    capacity = m->capacity ? 2 * m->capacity : 8;
    if (capacity > SIZE_MAX / sizeof(*rows)) return -1;
    rows = realloc(m->rows, capacity * sizeof(*rows));
    if (rows == NULL) return -1;
    m->rows = rows;
    m->capacity = capacity;
  }
  m->rows[m->nrows++] = *row;
  return 0;
}

void RestartManifestFree(struct RestartManifest *m)
{
  if (m == NULL) return;
  free(m->rows);
  memset(m, 0, sizeof(*m));
}

/* Pure generation selection (V1): no published filename may be overwritten before manifest publication. */
int RestartChooseGeneration(const struct RestartManifest *old, const uint64_t *sets, int nset, uint64_t *gen)
{
  for (uint64_t candidate = 0; candidate < 2; ++candidate) {
    int collision = 0;
    for (int i = 0; i < nset; ++i) {
      const struct RestartManifestRow *row = RestartManifestFind(old, sets[i]);
      if (row != NULL && row->gen == candidate) collision = 1;
    }
    if (!collision) { *gen = candidate; return 0; }
  }
  return -1;
}

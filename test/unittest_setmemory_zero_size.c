/* Unit test for zero-size 2d/3d setmemory allocators.
 * ASan catches the pre-fix A[0] / A[0][0] writes into calloc(0) storage. */
#include "setmemory.h"
#include <stdio.h>

static int g_failed = 0;

static void check_ptr(const char *name, const void *ptr) {
  if (ptr == NULL) {
    printf("%s : NULL\n", name);
    g_failed = 1;
  }
}

static void check_li_2d(unsigned long n, unsigned long m) {
  long int **a = li_2d_allocate(n, m);
  check_ptr("li_2d A", a);
  if (a != NULL) check_ptr("li_2d A[0]", a[0]);
  if (a != NULL && a[0] != NULL) free_li_2d_allocate(a);
}

static void check_i_2d(unsigned long n, unsigned long m) {
  int **a = i_2d_allocate(n, m);
  check_ptr("i_2d A", a);
  if (a != NULL) check_ptr("i_2d A[0]", a[0]);
  if (a != NULL && a[0] != NULL) free_i_2d_allocate(a);
}

static void check_d_2d(unsigned long n, unsigned long m) {
  double **a = d_2d_allocate(n, m);
  check_ptr("d_2d A", a);
  if (a != NULL) check_ptr("d_2d A[0]", a[0]);
  if (a != NULL && a[0] != NULL) free_d_2d_allocate(a);
}

static void check_cd_2d(unsigned long n, unsigned long m) {
  double complex **a = cd_2d_allocate(n, m);
  check_ptr("cd_2d A", a);
  if (a != NULL) check_ptr("cd_2d A[0]", a[0]);
  if (a != NULL && a[0] != NULL) free_cd_2d_allocate(a);
}

static void check_i_3d(unsigned long n, unsigned long m, unsigned long l) {
  int ***a = i_3d_allocate(n, m, l);
  check_ptr("i_3d A", a);
  if (a != NULL) check_ptr("i_3d A[0]", a[0]);
  if (a != NULL && a[0] != NULL) check_ptr("i_3d A[0][0]", a[0][0]);
  if (a != NULL && a[0] != NULL && a[0][0] != NULL) free_i_3d_allocate(a);
}

static void check_cd_3d(unsigned long n, unsigned long m, unsigned long l) {
  double complex ***a = cd_3d_allocate(n, m, l);
  check_ptr("cd_3d A", a);
  if (a != NULL) check_ptr("cd_3d A[0]", a[0]);
  if (a != NULL && a[0] != NULL) check_ptr("cd_3d A[0][0]", a[0][0]);
  if (a != NULL && a[0] != NULL && a[0][0] != NULL) free_cd_3d_allocate(a);
}

int main(void) {
  check_li_2d(0, 2);
  check_li_2d(2, 0);
  check_i_2d(0, 2);
  check_i_2d(2, 0);
  check_d_2d(0, 2);
  check_d_2d(2, 0);
  check_cd_2d(0, 2);
  check_cd_2d(2, 0);

  check_i_3d(0, 2, 3);
  check_i_3d(2, 0, 3);
  check_i_3d(2, 3, 0);
  check_cd_3d(0, 2, 3);
  check_cd_3d(2, 0, 3);
  check_cd_3d(2, 3, 0);

  printf("\n%s\n", g_failed ? "UNIT TEST FAILED" : "UNIT TEST PASSED");
  return g_failed;
}

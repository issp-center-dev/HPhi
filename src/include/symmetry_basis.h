#ifndef HPHI_SYMMETRY_BASIS_H
#define HPHI_SYMMETRY_BASIS_H

#include "Common.h"

struct BindStruct;
struct DefineList;
struct SymmetryMatvecPlan;

struct SymmetryBasisVector {
  unsigned long int rep_state;
  unsigned int orbit_size;
  unsigned int stabilizer_size;
  double norm;
  double complex stabilizer_character_sum;
  double diagonal;
};

struct SymmetryCanonicalResult {
  int found;
  unsigned long int basis_index;
  unsigned int op_rep_to_state;
  double complex phase;
};

struct SymmetryTransformResult {
  unsigned long int state;
  double complex amplitude;
};

struct SymmetryBasisRuntime {
  int enabled;
  unsigned int nsite;
  unsigned int group_order;
  unsigned long int full_dim;
  unsigned long int dim;
  unsigned long int capacity;
  struct SymmetryBasisVector *basis;
  unsigned long int rep_hash_size;
  unsigned long int *rep_hash_keys;
  unsigned long int *rep_hash_values;
  unsigned int *group_inverse;
  unsigned long int local_offset;
  unsigned long int local_dim;
  unsigned long long basis_raw_states;
  unsigned long long basis_representative_candidates;
  unsigned long long basis_compatible_survivors;
  unsigned long long basis_transform_calls;
  unsigned long long basis_orbit_metadata_calls;
  unsigned int basis_thread_count;
  unsigned long long basis_thread_raw_states_max;
  unsigned long long basis_thread_representative_candidates_max;
  unsigned long long basis_thread_compatible_survivors_max;
  unsigned long long basis_thread_transform_calls_max;
  unsigned long long basis_gather_entries;
  unsigned long long basis_gather_bytes;
  int *mpi_recvcounts;
  int *mpi_displs;
  double complex *mpi_full_v1;
  int matvec_mode;
  int vector_exchange_mode;
  struct SymmetryMatvecPlan *matvec_plan;
};

static inline const struct SymmetryBasisVector *SymmetryBasisLocalEntry(
    const struct SymmetryBasisRuntime *sym,
    unsigned long int local_index)
{
  unsigned long int global_index;
  if (sym == NULL || sym->enabled != TRUE || sym->basis == NULL ||
      local_index == 0UL || local_index > sym->local_dim ||
      sym->local_offset > sym->dim ||
      sym->local_dim > sym->dim - sym->local_offset) {
    return NULL;
  }
  global_index = sym->local_offset + local_index;
  if (global_index > sym->capacity) return NULL;
  return &sym->basis[global_index];
}

static inline const struct SymmetryBasisVector *
SymmetryBasisReplicatedGlobalEntry(
    const struct SymmetryBasisRuntime *sym,
    unsigned long int global_index)
{
  if (sym == NULL || sym->enabled != TRUE || sym->basis == NULL ||
      global_index == 0UL || global_index > sym->dim ||
      global_index > sym->capacity) {
    return NULL;
  }
  return &sym->basis[global_index];
}

int ValidateSymmetryGroupInput(const struct DefineList *def);
unsigned long int SymmetryApplyToSpinBits(unsigned long int state,
                                          const int *perm,
                                          unsigned int nsite);
int SymmetryApplyToState(const struct DefineList *def,
                         unsigned long int state,
                         unsigned int op,
                         struct SymmetryTransformResult *result);
int BuildSymmetryBasis(struct BindStruct *X);
int SymmetryCanonicalizeState(const struct BindStruct *X,
                              unsigned long int state,
                              struct SymmetryCanonicalResult *result);
int SymmetryCanonicalizeSpinState(const struct BindStruct *X,
                                  unsigned long int state,
                                  struct SymmetryCanonicalResult *result);
int ActivateSymmetryBasisDimension(struct BindStruct *X);
int SymmetryBasisGlobalToLocal(const struct SymmetryBasisRuntime *sym,
                               unsigned long int global_index,
                               unsigned long int *local_index);
int GetOwnedHamiltonianDiagonal(const struct BindStruct *X,
                                unsigned long int local_index,
                                double *diagonal);
int ValidateSymmetrySectorOptions(const struct BindStruct *X);
void FreeSymmetryBasis(struct SymmetryBasisRuntime *sym);

#endif /* HPHI_SYMMETRY_BASIS_H */

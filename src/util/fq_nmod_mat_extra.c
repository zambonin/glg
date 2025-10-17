#include "util/fq_nmod_mat_extra.h"

void fq_nmod_mat_randtest_not_zero(fq_nmod_mat_t v, const fq_nmod_ctx_t ctx,
                                   flint_rand_t state) {
  do {
    fq_nmod_mat_randtest(v, state, ctx);
  } while (fq_nmod_mat_is_zero(v, ctx));
}

slong fq_nmod_mat_first_pos_entry(fq_nmod_mat_t v, const fq_nmod_ctx_t ctx) {
  for (slong i = 0; i < v->r; ++i) {
    if (!fq_nmod_is_zero(fq_nmod_mat_entry(v, i, 0), ctx)) {
      return i;
    }
  }
  return -1;
}

void fq_nmod_mat_set_minor(fq_nmod_mat_t mat, const ulong skip_i,
                           const ulong skip_j, const fq_nmod_mat_t src,
                           const fq_nmod_ctx_t ctx) {
  const slong r = fq_nmod_mat_nrows(src, ctx);
  const slong c = fq_nmod_mat_ncols(src, ctx);

  for (slong i = 0; i < r; ++i) {
    for (slong j = 0; j < c; ++j) {
      fq_nmod_mat_entry_set(mat, i + (i >= skip_i), j + (j >= skip_j),
                            fq_nmod_mat_entry(src, i, j), ctx);
    }
  }
}

void fq_nmod_mat_scalar_mul(fq_nmod_mat_t B, const fq_nmod_mat_t A,
                            const fq_nmod_t c, const fq_nmod_ctx_t ctx) {
  for (slong i = 0; i < A->r; ++i) {
    for (slong j = 0; j < A->c; ++j) {
      fq_nmod_mul(fq_nmod_mat_entry(B, i, j), fq_nmod_mat_entry(A, i, j), c,
                  ctx);
    }
  }
}

void fmpz_randlimb_m(fmpz_t f, flint_rand_t state, const fmpz_t m) {
  fmpz_zero(f);
  const slong limbs = fmpz_size(m);
  ulong l[limbs];

  for (slong i = 0; i < limbs; ++i) {
    l[i] = n_randlimb(state);
  }

  fmpz_set_signed_ui_array(f, l, limbs);
  fmpz_mod(f, f, m);
}

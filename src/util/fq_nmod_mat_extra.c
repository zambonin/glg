#include "util/fq_nmod_mat_extra.h"
#include "comb/poly.h"
#include "util/rand.h"

void fq_nmod_mat_randtest_not_zero_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                                       bit_buffer_t *buf) {
  const slong r = M->r;
  const slong c = M->c;

  if (r == 0 || c == 0) {
    return;
  }

  fmpz_t q;
  fmpz_init(q);

  fmpz_t Q;
  fmpz_init(Q);

  fmpz_t sample;
  fmpz_init(sample);

  fmpz_t entry_rank;
  fmpz_init(entry_rank);

  fq_nmod_ctx_order(q, ctx);
  fmpz_pow_ui(Q, q, r * c);
  fmpz_sub_ui(Q, Q, 1);

  fast_dice_roller(sample, Q, buf);
  fmpz_add_ui(sample, sample, 1);

  for (slong i = 0; i < r; ++i) {
    for (slong j = 0; j < c; ++j) {
      fmpz_fdiv_qr(sample, entry_rank, sample, q);
      unrank_poly(fq_nmod_mat_entry(M, i, j), entry_rank, ctx);
    }
  }

  fmpz_clear(entry_rank);
  fmpz_clear(sample);
  fmpz_clear(Q);
  fmpz_clear(q);
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

void fq_nmod_mat_entry_rand_buf(fq_nmod_mat_t M, const slong i, const slong j,
                                const fq_nmod_ctx_t ctx, bit_buffer_t *buf) {
  fq_nmod_rand_buf(fq_nmod_mat_entry(M, i, j), ctx, buf);
}

void fq_nmod_mat_entry_rand_not_zero_buf(fq_nmod_mat_t M, const slong i,
                                         const slong j,
                                         const fq_nmod_ctx_t ctx,
                                         bit_buffer_t *buf) {
  fq_nmod_rand_not_zero_buf(fq_nmod_mat_entry(M, i, j), ctx, buf);
}

void fq_nmod_mat_randtest_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                              bit_buffer_t *buf) {
  for (slong i = 0; i < M->r; ++i) {
    for (slong j = 0; j < M->c; ++j) {
      fq_nmod_mat_entry_rand_buf(M, i, j, ctx, buf);
    }
  }
}

void fq_nmod_mat_randtriu_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                              bit_buffer_t *buf, const int unit) {
  const slong n = M->r;
  fq_nmod_mat_zero(M, ctx);

  for (slong i = 0; i < n; ++i) {
    for (slong j = i; j < n; ++j) {
      if (i == j) {
        if (unit) {
          fq_nmod_one(fq_nmod_mat_entry(M, i, j), ctx);
        } else {
          fq_nmod_mat_entry_rand_not_zero_buf(M, i, j, ctx, buf);
        }
      } else {
        fq_nmod_mat_entry_rand_buf(M, i, j, ctx, buf);
      }
    }
  }
}

void fq_nmod_mat_randtril_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                              bit_buffer_t *buf, const int unit) {
  fq_nmod_mat_randtriu_buf(M, ctx, buf, unit);
  fq_nmod_mat_transpose(M, M, ctx);
}

#include "comb/gl.h"
#include "comb/sspace.h"

void unrank_gl(fq_nmod_mat_t M, const fmpz_t r, const fq_nmod_ctx_t ctx) {
  const ulong n = fq_nmod_mat_nrows(M, ctx);
  fq_nmod_mat_zero(M, ctx);

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_n;
  fmpz_init(q_pow_n);
  fmpz_pow_ui(q_pow_n, q, n);

  fmpz_t q_pow_i;
  fmpz_init_set_ui(q_pow_i, 1);

  fmpz_t rp;
  fmpz_init_set(rp, r);

  fmpz_t count;
  fmpz_init(count);

  fmpz_t k;
  fmpz_init(k);

  fq_nmod_mat_t v;
  fq_nmod_mat_init(v, n, 1, ctx);

  for (ulong i = 0; i < n; ++i) {
    fmpz_zero(k);
    fmpz_sub(count, q_pow_n, q_pow_i);

    if (fmpz_cmp_ui(count, 0) > 0) {
      fmpz_fdiv_qr(rp, k, rp, count);
    } else {
      fmpz_set(k, rp);
    }

    fq_nmod_mat_t B;
    fq_nmod_mat_window_init(B, M, 0, 0, n, i, ctx);
    unrank_subspace(v, B, k, ctx);
    fq_nmod_mat_window_clear(B, ctx);

    for (ulong j = 0; j < n; ++j) {
      fq_nmod_set(fq_nmod_mat_entry(M, j, i), fq_nmod_mat_entry(v, j, 0), ctx);
    }

    fmpz_mul(q_pow_i, q_pow_i, q);
  }

  fq_nmod_mat_clear(v, ctx);
  fmpz_clear(k);
  fmpz_clear(count);
  fmpz_clear(rp);
  fmpz_clear(q_pow_i);
  fmpz_clear(q_pow_n);
  fmpz_clear(q);
}

void rank_gl(fmpz_t r, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx) {
  const ulong n = fq_nmod_mat_nrows(M, ctx);
  fmpz_zero(r);

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_n;
  fmpz_init(q_pow_n);
  fmpz_pow_ui(q_pow_n, q, n);

  fmpz_t q_pow_i;
  fmpz_init_set_ui(q_pow_i, 1);

  fmpz_t acc;
  fmpz_init_set_ui(acc, 1);

  fmpz_t count;
  fmpz_init(count);

  fmpz_t k;
  fmpz_init(k);

  fq_nmod_mat_t v;
  fq_nmod_mat_init(v, n, 1, ctx);

  for (ulong i = 0; i < n; ++i) {
    fmpz_zero(k);
    fmpz_zero(count);

    for (ulong j = 0; j < n; ++j) {
      fq_nmod_set(fq_nmod_mat_entry(v, j, 0), fq_nmod_mat_entry(M, j, i), ctx);
    }

    fq_nmod_mat_t B;
    fq_nmod_mat_window_init(B, M, 0, 0, n, i, ctx);
    rank_subspace(k, v, B, ctx);
    fq_nmod_mat_window_clear(B, ctx);

    fmpz_addmul(r, k, acc);
    fmpz_sub(count, q_pow_n, q_pow_i);
    fmpz_mul(acc, acc, count);
    fmpz_mul(q_pow_i, q_pow_i, q);
  }

  fq_nmod_mat_clear(v, ctx);
  fmpz_clear(k);
  fmpz_clear(count);
  fmpz_clear(acc);
  fmpz_clear(q_pow_i);
  fmpz_clear(q_pow_n);
  fmpz_clear(q);
}

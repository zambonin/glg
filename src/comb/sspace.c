#include "comb/sspace.h"
#include "comb/vspace.h"
#include "util/math.h"

void unrank_subspace(fq_nmod_mat_t v, const fq_nmod_mat_t B, const fmpz_t r,
                     const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(B, ctx);
  const slong n_i = fq_nmod_mat_ncols(B, ctx);
  const slong n_o = n - n_i;
  fq_nmod_mat_zero(v, ctx);

  if (n_o == 0) {
    return;
  }

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_n_o;
  fmpz_init(q_pow_n_o);
  fmpz_pow_ui(q_pow_n_o, q, n_o);
  fmpz_sub_ui(q_pow_n_o, q_pow_n_o, 1);

  fmpz_t r_i;
  fmpz_init(r_i);
  fmpz_t r_o;
  fmpz_init(r_o);
  if (fmpz_cmp_ui(q_pow_n_o, 0) > 0) {
    fmpz_fdiv_qr(r_i, r_o, r, q_pow_n_o);
  } else {
    fmpz_set(r_i, r);
    fmpz_zero(r_o);
  }

  fq_nmod_mat_t C_i;
  fq_nmod_mat_init(C_i, n_i, 1, ctx);

  fq_nmod_mat_t C_o;
  fq_nmod_mat_init(C_o, n_o, 1, ctx);

  fq_nmod_mat_t c;
  fq_nmod_mat_init(c, n, 1, ctx);

  fq_nmod_mat_t F;
  fq_nmod_mat_init(F, n, n, ctx);

  unrank_vspace(C_i, r_i, ctx);
  fmpz_add_ui(r_o, r_o, 1);
  unrank_vspace(C_o, r_o, ctx);

  complete_basis(F, B, n, ctx);

  for (ulong j = 0; j < n_i; ++j) {
    fq_nmod_set(fq_nmod_mat_entry(c, j, 0), fq_nmod_mat_entry(C_i, j, 0), ctx);
  }
  for (ulong j = n_i; j < n; ++j) {
    ulong l = j - n_i;
    fq_nmod_set(fq_nmod_mat_entry(c, j, 0), fq_nmod_mat_entry(C_o, l, 0), ctx);
  }

  fq_nmod_mat_mul(v, F, c, ctx);

  fq_nmod_mat_clear(F, ctx);
  fq_nmod_mat_clear(c, ctx);
  fq_nmod_mat_clear(C_o, ctx);
  fq_nmod_mat_clear(C_i, ctx);
  fmpz_clear(r_o);
  fmpz_clear(r_i);
  fmpz_clear(q_pow_n_o);
  fmpz_clear(q);
}

void rank_subspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_mat_t B,
                   const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(B, ctx);
  const slong n_i = fq_nmod_mat_ncols(B, ctx);
  const slong n_o = n - n_i;
  fmpz_zero(r);

  if (n_o == 0) {
    return;
  }

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_n_o;
  fmpz_init(q_pow_n_o);
  fmpz_pow_ui(q_pow_n_o, q, n_o);
  fmpz_sub_ui(q_pow_n_o, q_pow_n_o, 1);

  fmpz_t r_i;
  fmpz_init(r_i);

  fmpz_t r_o;
  fmpz_init(r_o);

  fq_nmod_mat_t C_i;
  fq_nmod_mat_init(C_i, n_i, 1, ctx);

  fq_nmod_mat_t C_o;
  fq_nmod_mat_init(C_o, n_o, 1, ctx);

  fq_nmod_mat_t c;
  fq_nmod_mat_init(c, n, 1, ctx);

  fq_nmod_mat_t F;
  fq_nmod_mat_init(F, n, n, ctx);

  complete_basis(F, B, n, ctx);
  (void)fq_nmod_mat_solve(c, F, v, ctx);

  for (ulong j = 0; j < n_i; ++j) {
    fq_nmod_set(fq_nmod_mat_entry(C_i, j, 0), fq_nmod_mat_entry(c, j, 0), ctx);
  }
  for (ulong j = n_i; j < n; ++j) {
    ulong l = j - n_i;
    fq_nmod_set(fq_nmod_mat_entry(C_o, l, 0), fq_nmod_mat_entry(c, j, 0), ctx);
  }

  rank_vspace(r_i, C_i, ctx);
  rank_vspace(r_o, C_o, ctx);
  fmpz_sub_ui(r_o, r_o, 1);

  if (fmpz_cmp_ui(q_pow_n_o, 0) > 0) {
    fmpz_mul(r, r_i, q_pow_n_o);
    fmpz_add(r, r, r_o);
  } else {
    fmpz_set(r, r_i);
  }

  fq_nmod_mat_clear(F, ctx);
  fq_nmod_mat_clear(c, ctx);
  fq_nmod_mat_clear(C_o, ctx);
  fq_nmod_mat_clear(C_i, ctx);
  fmpz_clear(r_o);
  fmpz_clear(r_i);
  fmpz_clear(q_pow_n_o);
  fmpz_clear(q);
}
